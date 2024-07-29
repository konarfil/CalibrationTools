#ifndef ENERGY_CORRECTION_CALCULATOR_H
#define ENERGY_CORRECTION_CALCULATOR_H

// Other
#include <vector>
#include <string>

class energy_correction_calculator
{
  public:
    enum material
    {
      gas,
      mylar,
      nylon
    };

    // Constructor
    energy_correction_calculator(double gas_pressure, double He_pressure, double Et_pressure,
                                 double Ar_pressure, double T_gas, std::string uniformity_correction_parameters_mwall_8inch_path,
                                 std::string uniformity_correction_parameters_mwall_5inch_path,
                                 std::string uniformity_correction_parameters_xwall_path);
  
    // Destructor
    ~energy_correction_calculator();

    // calculate energy losses in material using the Bethe-Bloch formula for electrons
    double Bethe_Bloch_loss(double d, double E, material mat);

    // based on observed energy and geometrical correction factor 
    // calculates energy corrected by the non-linearity effects
    double non_linearity_correction(double Ef, double geometrical_factor);

    // calculate geometrical correction factor based on vertex position relative to the OM and OM_type
    double geometrical_factor(double x, double y, double z, int OM_type);

  private:

    // parse parameters for the non-uniformity correction
    std::vector<double> parse_pol3d_parameters(const std::string & parameters_path);

    // intensities of Bi207 internal conversion electrons
    static constexpr double I_482keV_  = 1.52;
    static constexpr double I_555keV_  = 0.44;
    static constexpr double I_567keV_  = 0.15;
    static constexpr double I_976keV_  = 7.03;
    static constexpr double I_1049keV_ = 1.84;
    static constexpr double I_1061keV_ = 0.54;

    static constexpr double E_0_ = 511.0; // rest energy of electron

    // properties of Mylar and nylon for energy loss calculation using the Bethe-Bloch formula
    static constexpr double Z_A_mylar_ = 0.520379; // Z/A ratio
    static constexpr double I_mylar_   = 0.0787; // mean ionization energy keV
    static constexpr double rho_mylar_ = 1.4; // density g/cm^3
    static constexpr double Z_A_nylon_ = 0.547912; // Z/A ratio
    static constexpr double I_nylon_   = 0.0639; // mean ionization energy keV
    static constexpr double rho_nylon_ = 1.14; // density g/cm^3

    // atomic (molecular) masses of gas components
    static constexpr double Ar_He_ = 4.002602;
    static constexpr double Ar_Et_ = 46.068;
    static constexpr double Ar_H_  = 1.00784;
    static constexpr double Ar_C_  = 12.011;
    static constexpr double Ar_O_  = 15.999;
    static constexpr double Ar_Ar_ = 39.948;

    // atomic numbers of gas components
    static const int Z_He_ = 2;
    static const int Z_H_  = 1;
    static const int Z_C_  = 6;
    static const int Z_O_  = 8;
    static const int Z_Ar_ = 18;

    // mean ionization energies of gas components in keV
    static constexpr double I_He_ = 0.0418;
    static constexpr double I_H_  = 0.0192;
    static constexpr double I_C_  = 0.0780;
    static constexpr double I_O_  = 0.095;
    static constexpr double I_Ar_ = 0.188;

    // gas properites for energy loss calculation using the Bethe-Bloch formula
    double Z_A_gas_;
    double I_gas_;
    double rho_gas_;

    std::vector<double> uniformity_correction_parameters_mwall_8inch_; // Polynomial parameters for the uniformity correction for MWall 8"
    std::vector<double> uniformity_correction_parameters_mwall_5inch_; // Polynomial parameters for the uniformity correction for MWall 5"
    std::vector<double> uniformity_correction_parameters_xwall_;       // Polynomial parameters for the uniformity correction for XWall

};

#endif