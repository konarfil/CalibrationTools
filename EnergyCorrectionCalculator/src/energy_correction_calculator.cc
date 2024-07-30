#include "energy_correction_calculator.h"

// Falaise
#include <falaise/snemo/processing/mock_calorimeter_s2c_module_utils.h>

// Other
#include <cmath>
#include <fstream>
#include <numeric>

// Constructor
energy_correction_calculator::energy_correction_calculator(double gas_pressure, double He_pressure, double Et_pressure,
                                 double Ar_pressure, double T_gas, std::string uniformity_correction_parameters_mwall_8inch_path,
                                 std::string uniformity_correction_parameters_mwall_5inch_path,
                                 std::string uniformity_correction_parameters_xwall_path)
{
  // calculate mass fractions of gas components
  double frac_He = (He_pressure * Ar_He_) / (He_pressure * Ar_He_ + Et_pressure * Ar_Et_ + Ar_pressure * Ar_Ar_);
  double frac_Et = (Et_pressure * Ar_Et_) / (He_pressure * Ar_He_ + Et_pressure * Ar_Et_ + Ar_pressure * Ar_Ar_);
  double frac_Ar = (Ar_pressure * Ar_Ar_) / (He_pressure * Ar_He_ + Et_pressure * Ar_Et_ + Ar_pressure * Ar_Ar_);
  double frac_C  = (frac_Et * 2 * Ar_C_) / (2*Ar_C_ + 6*Ar_H_ + Ar_O_);
  double frac_O  = (frac_Et * Ar_O_) / (2*Ar_C_ + 6*Ar_H_ + Ar_O_);
  double frac_H  = (frac_Et * 6 * Ar_H_) / (2*Ar_C_ + 6*Ar_H_ + Ar_O_);

  // calculate Z/A ratio of the tracking gas
  Z_A_gas_ = (frac_He * Z_He_) / Ar_He_ + (frac_H * Z_H_) / Ar_H_ 
            + (frac_C * Z_C_) / Ar_C_ + (frac_O * Z_O_) / Ar_O_ + (frac_Ar * Z_Ar_) / Ar_Ar_;

  // calculate mean ionization energy of the tracking gas
  I_gas_ = std::exp(((frac_He * Z_He_ * std::log(I_He_)) / Ar_He_ 
                   + (frac_H * Z_H_ * std::log(I_H_)) / Ar_H_
                   + (frac_C * Z_C_ * std::log(I_C_)) / Ar_C_ 
                   + (frac_O * Z_O_ * std::log(I_O_)) / Ar_O_ 
                   + (frac_Ar * Z_Ar_ * std::log(I_Ar_)) / Ar_Ar_) / Z_A_gas_);

  // calculate molar mas of the tracking gas
  double Mm_gas = He_pressure * Ar_He_ + Et_pressure * Ar_Et_ + Ar_pressure * Ar_Ar_;

  // calculate density of the tracking gas (g/cm^3)
  rho_gas_ = Mm_gas / (8.31 * T_gas) * gas_pressure * 0.1;

  // Initialize pol3d parameters for non-uniformity correction
  uniformity_correction_parameters_mwall_8inch_ = parse_pol3d_parameters(uniformity_correction_parameters_mwall_8inch_path);
  uniformity_correction_parameters_mwall_5inch_ = parse_pol3d_parameters(uniformity_correction_parameters_mwall_5inch_path);
  uniformity_correction_parameters_xwall_ = parse_pol3d_parameters(uniformity_correction_parameters_xwall_path);
}

// Destructor
energy_correction_calculator::~energy_correction_calculator()
{

}

// parse parameters for the non-uniformity correction
std::vector<double> energy_correction_calculator::parse_pol3d_parameters(const std::string & parameters_path)
{
  std::ifstream parameters_file(parameters_path.c_str());

  double par, par_err;
  std::vector<double> parameters;

  while(parameters_file >> par >> par_err)
    parameters.push_back(par);

  return parameters;
}

// calculate energy losses in material using the Bethe-Bloch formula for electrons
double energy_correction_calculator::Bethe_Bloch_loss(double d, double E, material mat)
{
  double Z_A, I, RHO;
  switch(mat)
  {
    case material::gas:
      Z_A = Z_A_gas_;
      I = I_gas_;
      RHO = rho_gas_;
      break;
    case material::mylar:
      Z_A = Z_A_mylar_;
      I = I_mylar_;
      RHO = rho_mylar_;
      break;
    case material::nylon:
      Z_A = Z_A_nylon_;
      I = I_nylon_;
      RHO = rho_nylon_;
      break;
  }

  double C = 307.075 / 2.0;
  
  double beta2 = (2.0 * E_0_ * E + E * E) / ((E_0_ + E) * (E_0_ + E));
  double stopping_power = C * Z_A / beta2 * (std::log((E_0_ * beta2 * E) / (2.0 * I * I * (1.0 - beta2))) 
  			- (2.0 * sqrt(1.0 - beta2) - 1.0 + beta2) * std::log(2.0) + 1.0 - beta2 + 1.0 / 8.0 * (1.0 - sqrt(1.0 - beta2)) * (1.0 - sqrt(1.0 - beta2)));
  double delta_E = stopping_power * RHO * d * 1e-1;
  
  return delta_E;
}

// based on observed energy and geometrical correction factor 
// calculates energy corrected by the non-linearity effects
// uses Newton's method to find the energy numerically
double energy_correction_calculator::non_linearity_correction(double Ef, double geometrical_factor)
{
  double Ei_prev = 1e10;
  double Ei = Ef;
  int i = 0;
  while(std::abs(Ei - Ei_prev) > 0.001 && i < 20)
  {
  	double f = 1.092096 * Ei - 1.56416*std::pow(Ei, 0.59) - Ef / geometrical_factor;
  	double df = 1.097086 - 0.9228544/std::pow(Ei, 0.41);
  	Ei_prev = Ei;
  	Ei = Ei - f/df;
  	i++;
  }
  return Ei;
}

// calculate geometrical correction factor based on vertex position relative to the OM and OM_type
double energy_correction_calculator::geometrical_factor(double x, double y, double z, int OM_type)
{
  double position_xyz[3] = {x, y, z};
  switch (OM_type)
  {
    case 1302: // M-wall
      position_xyz[2] += 15.50000001; // add half height of scintillator
      if (std::abs(z) < 1500.)
        return snemo::processing::pol3d(position_xyz, &uniformity_correction_parameters_mwall_8inch_[0]);
      else
        return snemo::processing::pol3d(position_xyz, &uniformity_correction_parameters_mwall_5inch_[0]);
      break;

    case 1232: // X-wall
      position_xyz[2] += 75.10000001; // add half height of scintillator
      return snemo::processing::pol3d(position_xyz, &uniformity_correction_parameters_xwall_[0]);
      break;
  }
}