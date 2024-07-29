#include "calibration_parameter_finder.h"

// Falaise
#include <falaise/snemo/processing/mock_calorimeter_s2c_module_utils.h>

// ROOT
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>

// Other
#include <cmath>
#include <fstream>
#include <numeric>

// Constructor
calibration_parameter_finder::calibration_parameter_finder(TTree* OM_data, double gas_pressure, double He_pressure, double Et_pressure, double Ar_pressure, 
                                                           double T_gas, double minimization_threshold, int max_iterations, bool save_fitted_spectra)
{
  minimization_threshold_ = minimization_threshold;
  max_iterations_ = max_iterations;
  save_fitted_spectra_ = save_fitted_spectra;
  fitted_spectra_dir_ = FITTED_SPECTRA_DIR;
  best_calib_ = calib_info{0.0, 0.0, 0.0, 0.0, 1e9};

  // set branch addresses
  OM_data_ = OM_data;
  source_vertex_pos_ = 0;
  calo_vertex_pos_ = 0;
  calo_vertex_pos_OM_ = 0;
  OM_data_->SetBranchAddress("charge", &charge_);
  OM_data_->SetBranchAddress("source_vertex_pos", &source_vertex_pos_);
  OM_data_->SetBranchAddress("calo_vertex_pos", &calo_vertex_pos_);
  OM_data_->SetBranchAddress("calo_vertex_pos_OM", &calo_vertex_pos_OM_);

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

  std::string falaise_resources_path = Falaise_RESOURCE_DIR;
  // Initialize the pol3d parameters for MWall 8"
  std::string pol3d_parameters_mwall_8inch_path = falaise_resources_path + "/snemo/demonstrator/reconstruction/db/fit_parameters_10D_MW_8inch.db";
  uniformity_correction_parameters_mwall_8inch_ = this->parse_pol3d_parameters(pol3d_parameters_mwall_8inch_path);

  // Initialize the pol3d parameters for MWall 5"
  std::string pol3d_parameters_mwall_5inch_path = falaise_resources_path + "/snemo/demonstrator/reconstruction/db/fit_parameters_10D_MW_5inch.db";
  uniformity_correction_parameters_mwall_5inch_ = this->parse_pol3d_parameters(pol3d_parameters_mwall_5inch_path);

  // Initialize the pol3d parameters for XWall
  std::string pol3d_parameters_xwall_path = falaise_resources_path + "/snemo/demonstrator/reconstruction/db/fit_parameters_10D_XW.db";
  uniformity_correction_parameters_xwall_ = this->parse_pol3d_parameters(pol3d_parameters_xwall_path);
}

// Destructor
calibration_parameter_finder::~calibration_parameter_finder()
{

}

// calculate calibration parameters of the OM and assign them to a and b
// the method minimizes the loss function using the Downhill Simplex Method
calib_info calibration_parameter_finder::find_calibration_parameters()
{
  // set initial points of the simplex
  double simplex[N_minimization_points_][N_calibration_parameters_] = {{100.0, -50.0},
                                                                      {150.0, 50.0},
                                                                      {200.0, -50.0}};

  // calculate loss function for the simplex points
  double loss[N_minimization_points_];
  for(int i = 0; i < N_minimization_points_; i++)
  {
    loss[i] = loss_function(simplex[i][0], simplex[i][1]);
  }
  
  for(int i = 0; i < max_iterations_; i++)
  {
    // Order the simplex vertices
    int indices[N_minimization_points_];
    std::iota(indices, indices + N_minimization_points_, 0);
    std::sort(&indices[0], &indices[N_minimization_points_], [&](int a, int b) { return loss[a] < loss[b]; });

    int best = indices[0];
    int worst = indices[N_minimization_points_ - 1];
    int secondWorst = indices[N_minimization_points_ - 2];

    // Check for convergence
    if (loss[best] < minimization_threshold_ || i == max_iterations_ - 1)
    {
      best_calib_.a = simplex[best][0];
      best_calib_.b = simplex[best][1];

      if(save_fitted_spectra_) // save pngs of the fitted spectra
      {
        TCanvas* c = new TCanvas("c", "c", 1280, 720);
        min_loss_histo_.GetXaxis()->SetRangeUser(0.0, 2000.0);
        min_loss_histo_.Draw();
        min_loss_gaussA_.Draw("same");
        min_loss_gaussB_.Draw("same");
        std::string gid_str = OM_data_->GetName();
        c->SaveAs((fitted_spectra_dir_ + "/" + gid_str + ".png").c_str());
        delete c;
      }
      return best_calib_;
    }

    // Calculate the centroid of the best n points
    double c[N_calibration_parameters_] = {(simplex[best][0] + simplex[secondWorst][0]) / 2.0,
                                           (simplex[best][1] + simplex[secondWorst][1]) / 2.0};

    // Reflection
    double xr[N_calibration_parameters_];
    for (int j = 0; j < N_calibration_parameters_; ++j)
    {
        xr[j] = c[j] + ALPHA_ * (c[j] - simplex[worst][j]);
    }
    double fxr = loss_function(xr[0], xr[1]);
    if (loss[best] <= fxr && fxr < loss[secondWorst])
    {
        std::copy(xr, xr + N_calibration_parameters_, simplex[worst]);
        loss[worst] = fxr;
        continue;
    }

    // Expansion
    if (fxr < loss[best])
    {
        double xe[N_calibration_parameters_];
        for (int j = 0; j < N_calibration_parameters_; ++j)
        {
            xe[j] = c[j] + GAMMA_ * (xr[j] - c[j]);
        }
        double fxe = loss_function(xe[0], xe[1]);
        if (fxe < fxr)
        {
            std::copy(xe, xe + N_calibration_parameters_, simplex[worst]);
            loss[worst] = fxe;
        } 
        else
        {
            std::copy(xr, xr + N_calibration_parameters_, simplex[worst]);
            loss[worst] = fxr;
        }
        continue;
    }

    // Contraction
    double xc[N_calibration_parameters_];
    for (int j = 0; j < N_calibration_parameters_; ++j)
    {
        xc[j] = c[j] + RHO_ * (simplex[worst][j] - c[j]);
    }
    double fxc = loss_function(xc[0], xc[1]);
    if (fxc < loss[worst])
    {
        std::copy(xc, xc + N_calibration_parameters_, simplex[worst]);
        loss[worst] = fxc;
        continue;
    }

    // Shrink
    for (int i = 1; i < N_minimization_points_; ++i) // Skip the best point
    {
      for (int j = 0; j < N_calibration_parameters_; ++j)
      {
          simplex[i][j] = simplex[best][j] + SIGMA_ * (simplex[i][j] - simplex[best][j]);
      }
      loss[i] = loss_function(simplex[i][0], simplex[i][1]);
    }
  }
}

// fitting function for the first Bi207 peak
double calibration_parameter_finder::tripleGaussA(double* x,  double* par)
{
  double PDF    = 0.0;
  double gauss1 = 0.0;
  double gauss2 = 0.0;
  double gauss3 = 0.0;
  
  // fitting parameters
  double n1   = par[0];
  double mu1  = par[1];
  double sig1 = par[2];
  double C    = par[3];
  
  // fix parameters of second and third gaussian based intensities and energies of Bi207 internal conversion electrons
  double n2   = (I_555keV_ / I_482keV_) * n1;
  double mu2  = (555.0 / 482.0) * mu1;
  double sig2 = sqrt(555.0 / 482.0) * sig1;
  double n3   = (I_567keV_ / I_482keV_) * n1;
  double mu3  = (567.0 / 482.0) * mu1;
  double sig3 = sqrt(567.0 / 482.0) * sig1;
  
  // Calculation of exponents of each Gaussian
  double arg1 = (sig1 != 0.0) ? (x[0] - mu1) / (sig1) : 0.0;
  double arg2 = (sig2 != 0.0) ? (x[0] - mu2) / (sig2) : 0.0;
  double arg3 = (sig3 != 0.0) ? (x[0] - mu3) / (sig3) : 0.0;
  
  // Calculate gaussians
  gauss1 = exp(-0.5 * arg1 * arg1) / ( sig1 * sqrt(2.0 * TMath::Pi()));
  gauss2 = exp(-0.5 * arg2 * arg2) / ( sig2 * sqrt(2.0 * TMath::Pi()));
  gauss3 = exp(-0.5 * arg3 * arg3) / ( sig3 * sqrt(2.0 * TMath::Pi()));
  
  // we add the C parameter as the first peak seems to lie on the "tail" of the second peak
  // maybe could be replaced by different model of the spectrum, but this seems to be working quite well
  PDF = n1*gauss1 + n2*gauss2 + n3*gauss3 + C;
  
  return PDF;
}

// fitting function for the second Bi207 peak
double calibration_parameter_finder::tripleGaussB(double* x,  double* par)
{
  double PDF    = 0.0;
  double gauss1 = 0.0;
  double gauss2 = 0.0;
  double gauss3 = 0.0;
  
  // fitting parameters
  double n1   = par[0];
  double mu1  = par[1];
  double sig1 = par[2];
  
  // fix parameters of second and third gaussian based intensities and energies of Bi207 internal conversion electrons
  double n2   = (I_1049keV_ / I_976keV_) * n1;
  double mu2  = (1049.0 / 976.0) * mu1;
  double sig2 = sqrt(1049.0 / 976.0) * sig1;
  double n3   = (I_1061keV_ / I_976keV_) * n1;
  double mu3  = (1061.0 / 976.0) * mu1;
  double sig3 = sqrt(1061.0 / 976.0) * sig1;
  
  // Calculation of exponents of each Gaussian
  double arg1 = (sig1 != 0.0) ? (x[0] - mu1) / (sig1) : 0.0;
  double arg2 = (sig2 != 0.0) ? (x[0] - mu2) / (sig2) : 0.0;
  double arg3 = (sig3 != 0.0) ? (x[0] - mu3) / (sig3) : 0.0;
  
  // Calculate gaussians
  gauss1 = exp(-0.5 * arg1 * arg1) / ( sig1 * sqrt(2.0 * TMath::Pi()));
  gauss2 = exp(-0.5 * arg2 * arg2) / ( sig2 * sqrt(2.0 * TMath::Pi()));
  gauss3 = exp(-0.5 * arg3 * arg3) / ( sig3 * sqrt(2.0 * TMath::Pi()));
  
  PDF = n1*gauss1 + n2*gauss2 + n3*gauss3;
  
  return PDF;
}

// parse parameters for the non-uniformity correction
std::vector<double> calibration_parameter_finder::parse_pol3d_parameters(const std::string & parameters_path)
{
  std::ifstream parameters_file(parameters_path.c_str());

  double par, par_err;
  std::vector<double> parameters;

  while(parameters_file >> par >> par_err)
    parameters.push_back(par);

  return parameters;
}

// calculate energy losses in material using the Bethe-Bloch formula for electrons
double calibration_parameter_finder::Bethe_Bloch_loss(double d, double E, double Z_A, double I, double RHO)
{
  double C = 307.075 / 2.0;
  
  double beta2 = (2.0 * E_0_ * E + E * E) / ((E_0_ + E) * (E_0_ + E));
  double stopping_power = C * Z_A / beta2 * (TMath::Log((E_0_ * beta2 * E) / (2.0 * I * I * (1.0 - beta2))) 
  			- (2.0 * sqrt(1.0 - beta2) - 1.0 + beta2) * TMath::Log(2.0) + 1.0 - beta2 + 1.0 / 8.0 * (1.0 - sqrt(1.0 - beta2)) * (1.0 - sqrt(1.0 - beta2)));
  double delta_E = stopping_power * RHO * d * 1e-1;
  
  return delta_E;
}

// based on observed energy and geometrical correction factor 
// calculates energy corrected by the non-linearity effects
double calibration_parameter_finder::non_linearity_correction(double Ef, double geometrical_factor)
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

double calibration_parameter_finder::loss_function(double a, double b)
{
  TH1F histo = TH1F("histo", "histo", 2000, -8000.0, 12000.0);
  const int N_entries = OM_data_->GetEntries();
  

  //iterate through all electrons, apply corrections and compose energy spectrum
  for(int i = 0;i < N_entries;i++)
  {
    OM_data_->GetEntry(i);

    // observed energy
  	double Ef = a * charge_ + b;

    // determine the factor for geometrical non-uniformity correction
    double position_xyz[3] = {calo_vertex_pos_OM_->X(), calo_vertex_pos_OM_->Y(), calo_vertex_pos_OM_->Z()};
	  double geometrical_factor = 1.0;
    std::string gid_str = OM_data_->GetName();
    int OM_type = std::stoi(gid_str.substr(0, 4));
	  switch (OM_type)
	  {
      case 1302: // M-wall
        position_xyz[2] += 15.50000001; // add half height of scintillator
        if (std::abs(calo_vertex_pos_OM_->Z()) < 1500.)
          geometrical_factor = snemo::processing::pol3d(position_xyz, &uniformity_correction_parameters_mwall_8inch_[0]);
        else
          geometrical_factor = snemo::processing::pol3d(position_xyz, &uniformity_correction_parameters_mwall_5inch_[0]);
        break;

      case 1232: // X-wall
        position_xyz[2] += 75.10000001; // add half height of scintillator
        geometrical_factor = snemo::processing::pol3d(position_xyz, &uniformity_correction_parameters_xwall_[0]);
        break;
	  }

    // distances passed through gas, Mylar and nylon
  	double dist_gas   = std::sqrt((calo_vertex_pos_->X() - source_vertex_pos_->X())*(calo_vertex_pos_->X() - source_vertex_pos_->X())
                                + (calo_vertex_pos_->Y() - source_vertex_pos_->Y())*(calo_vertex_pos_->Y() - source_vertex_pos_->Y())
                                + (calo_vertex_pos_->Z() - source_vertex_pos_->Z())*(calo_vertex_pos_->Z() - source_vertex_pos_->Z()));

    double Ef_optical = non_linearity_correction(Ef, geometrical_factor); // apply optical correction
    double E;
    if(OM_type == 1302) // main wall
    {
      // 1 devided by the tracke angle measured relatively to the source foil normal
      double one_over_cos = dist_gas / std::abs(source_vertex_pos_->X() - calo_vertex_pos_->X());

      // calculate distance passed in mylar and nylon
      double dist_mylar = one_over_cos * mylar_thickness_;
      double dist_nylon = one_over_cos * nylon_thickness_;

      // apply energy loss correction
      double E_mylar = Ef_optical + Bethe_Bloch_loss(dist_mylar, Ef_optical, Z_A_mylar_, I_mylar_, rho_mylar_);
    	double E_nylon = E_mylar + Bethe_Bloch_loss(dist_nylon, E_mylar, Z_A_nylon_, I_nylon_, rho_nylon_);
    	E = E_nylon + Bethe_Bloch_loss(dist_gas, E_nylon, Z_A_gas_, I_gas_, rho_gas_);
    }
    else // xcalo
    {
      // 1 devided by the tracke angle measured relatively to the source foil normal
      double one_over_cos = dist_gas / std::abs(source_vertex_pos_->Y() - calo_vertex_pos_->Y());

      // calculate distance passed in mylar and nylon
      double dist_mylar = one_over_cos * mylar_thickness_;

      // apply energy loss correction (there is no nylon on xwall)
      double E_mylar = Ef_optical + Bethe_Bloch_loss(dist_mylar, Ef_optical, Z_A_mylar_, I_mylar_, rho_mylar_);
    	E = E_mylar + Bethe_Bloch_loss(dist_gas, E_mylar, Z_A_gas_, I_gas_, rho_gas_);
    }

  	histo.Fill(E);
  }

  //find two local maxima coresponding to 482 and 976 keV
	double histo_max   = histo.GetMaximum();
	double bin_width   = histo.GetBinWidth(1);
	double X_min       = histo.GetXaxis()->GetXmin();
	int    N_bins      = histo.GetNbinsX();
	
  // we find the fitting ranges using running average of the spectrum
	TH1F running_average = TH1F("avg", "avg", N_bins, histo.GetXaxis()->GetXmin(), histo.GetXaxis()->GetXmax());
	for(int i = 3;i < N_bins - 3;i++)
	{
		double sum = 0.0;
		for(int j = 0;j < 5;j++)
		{
			sum += histo.GetBinContent(i + j - 2);
		}
		running_average.SetBinContent(i, sum / 5.0);
	}
	
	double max_bin2 = running_average.GetMaximumBin();
	double max2     = running_average.GetBinContent(max_bin2);
	
	running_average.GetXaxis()->SetRangeUser(X_min, (X_min + max_bin2 * bin_width) * (2.0 / 3.0));
	double max_bin1 = running_average.GetMaximumBin();
	double max1     = running_average.GetBinContent(max_bin1);
	running_average.GetXaxis()->UnZoom();
	
  // find fitting ranges
  // the chosen conditions for fitting ranges are rather arbitrary but tested to work well
	double ranges[4] = {0.0, 0.0, 0.0, 0.0};
	for(int i = max_bin1;i > 1;i--)
	{			
		if(running_average.GetBinContent(i) < 0.5 * max1)
		{
			ranges[0] = X_min + i * bin_width;
			break;
		}
	}
	
	for(int i = max_bin1;i < N_bins;i++)
	{
		if(running_average.GetBinContent(i) < 0.5 * max1 || X_min + i * bin_width > 7.0/12.0*(X_min + max_bin2 * bin_width))
		{
			ranges[1] = X_min + i * bin_width;
			break;
		}
	}
	
	for(int i = max_bin2;i > 1;i--)
	{
		if(running_average.GetBinContent(i) < 0.20 * max2)
		{
			ranges[2] = X_min + i * bin_width;
			break;
		}
	}
	
	for(int i = max_bin2;i < N_bins;i++)
	{
		if(running_average.GetBinContent(i) < 0.15 * max2)
		{
			ranges[3] = X_min + i * bin_width;
			break;
		}
	}
	
  //define fitting functions; values of fitting parameters were chosen by trial and error
	TF1 tripGFitA = TF1("TripleGA", tripleGaussA, ranges[0], ranges[1], 4);
	TF1 tripGFitB = TF1("TripleGB", tripleGaussB, ranges[2], ranges[3], 3);
	tripGFitA.SetNpx(10000);
	tripGFitB.SetNpx(10000);
	tripGFitA.SetParNames("N1_A",  "mu1_A", "sig1_A", "C");
	tripGFitB.SetParNames("N1_B",  "mu1_B", "sig1_B");
	double max_bin1_x = histo.GetXaxis()->GetBinCenter(max_bin1);
	double sigma_est  = (max_bin1_x - ranges[0]) / 2.355;
	tripGFitA.SetParameters(histo_max * 75.0, max_bin1_x, sigma_est, 1.0);
	tripGFitA.SetParLimits(0, 0.0, 999999999.0);
	tripGFitA.SetParLimits(1, 0.9*max_bin1_x, 1.1*max_bin1_x);
	tripGFitA.SetParLimits(2, 0.0, 2.0*sigma_est);
	tripGFitA.SetParLimits(3, 0.0, 100.0);
	tripGFitB.SetParameters(histo_max * 360.0, (ranges[2] + ranges[3]) / 2.0, 0.8 * (ranges[3] - ranges[2]) / 2.0);
	tripGFitA.SetParLimits(0, 0.0, 999999999.0);
	
	histo.Fit(&tripGFitA, "QR");
	histo.Fit(&tripGFitB, "QR");
	
  //attempt multiple fits when chi-square is bad
	int i = 0;
	while(tripGFitA.GetChisquare() / tripGFitA.GetNDF() > 3.0 && i < 5)
	{
		histo.Fit(&tripGFitA, "QR");
		i++;
	}
	i = 0;
	while(tripGFitB.GetChisquare() / tripGFitB.GetNDF() > 3.0 && i < 5)
	{
		histo.Fit(&tripGFitB, "QR");
		i++;
	}
	
	double loss = (tripGFitA.GetParameter(1) - 482.0)*(tripGFitA.GetParameter(1) - 482.0) + (tripGFitB.GetParameter(1) - 976.0)*(tripGFitB.GetParameter(1) - 976.0);

  if(loss < best_calib_.loss)
  {
    best_calib_.loss = loss;
    best_calib_.chi2_NDF_A = tripGFitA.GetChisquare() / tripGFitA.GetNDF();
    best_calib_.chi2_NDF_B = tripGFitB.GetChisquare() / tripGFitB.GetNDF();
    if(save_fitted_spectra_)
    {
      min_loss_histo_ = histo;
      min_loss_histo_.SetName("min_loss_histo");
      min_loss_gaussA_ = tripGFitA;
      min_loss_gaussB_ = tripGFitB;
    }
  }

	return loss;
}