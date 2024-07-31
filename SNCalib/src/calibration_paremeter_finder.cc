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
calibration_parameter_finder::calibration_parameter_finder(TTree* OM_data, double minimization_threshold, int max_iterations,
                                                           bool save_fitted_spectra, energy_correction_calculator* corr_calculator)
{
  minimization_threshold_ = minimization_threshold;
  max_iterations_ = max_iterations;
  save_fitted_spectra_ = save_fitted_spectra;
  fitted_spectra_dir_ = FITTED_SPECTRA_DIR;
  best_calib_ = calib_info{0.0, 0.0, 0.0, 0.0, 1e9};

  corr_calculator_ = corr_calculator;

  // set branch addresses
  OM_data_ = OM_data;
  source_vertex_pos_ = 0;
  calo_vertex_pos_ = 0;
  calo_vertex_pos_OM_ = 0;
  OM_data_->SetBranchAddress("charge", &charge_);
  OM_data_->SetBranchAddress("source_vertex_pos", &source_vertex_pos_);
  OM_data_->SetBranchAddress("calo_vertex_pos", &calo_vertex_pos_);
  OM_data_->SetBranchAddress("calo_vertex_pos_OM", &calo_vertex_pos_OM_);
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

double calibration_parameter_finder::loss_function(double a, double b)
{
  TH1F histo = TH1F("histo", "histo", 2000, -8000.0, 12000.0);
  const int N_entries = OM_data_->GetEntries();
  

  //iterate through all electrons, apply corrections and compose energy spectrum
  for(int i = 0;i < N_entries;i++)
  {
    OM_data_->GetEntry(i);

    std::string om_num_str = OM_data_->GetName();
    int om_num = std::stoi(om_num_str);
    int OM_type;

    if(om_num < 520) // mwall
      OM_type = 1302;
    else // xcalo
      OM_type = 1232;

    // observed energy
  	double Ef = a * charge_ + b;

    // distance passed in gas
  	double dist_gas   = std::sqrt((calo_vertex_pos_->X() - source_vertex_pos_->X())*(calo_vertex_pos_->X() - source_vertex_pos_->X())
                                + (calo_vertex_pos_->Y() - source_vertex_pos_->Y())*(calo_vertex_pos_->Y() - source_vertex_pos_->Y())
                                + (calo_vertex_pos_->Z() - source_vertex_pos_->Z())*(calo_vertex_pos_->Z() - source_vertex_pos_->Z()));

    // apply optical correction
    double geometrical_factor = corr_calculator_->geometrical_factor(calo_vertex_pos_OM_->X(), calo_vertex_pos_OM_->Y(), calo_vertex_pos_OM_->Z(), OM_type);
    double Ef_optical = corr_calculator_->non_linearity_correction(Ef, geometrical_factor);
    double E;
    if(OM_type == 1302) // main wall
    {
      // 1 devided by the tracke angle measured relatively to the source foil normal
      double one_over_cos = dist_gas / std::abs(source_vertex_pos_->X() - calo_vertex_pos_->X());

      // calculate distance passed in mylar and nylon
      double dist_mylar = one_over_cos * mylar_thickness_;
      double dist_nylon = one_over_cos * nylon_thickness_;

      // apply energy loss correction
      double E_mylar = Ef_optical + corr_calculator_->Bethe_Bloch_loss(dist_mylar, Ef_optical, energy_correction_calculator::material::mylar);
    	double E_nylon = E_mylar + corr_calculator_->Bethe_Bloch_loss(dist_nylon, E_mylar, energy_correction_calculator::material::nylon);
    	E = E_nylon + corr_calculator_->Bethe_Bloch_loss(dist_gas, E_nylon, energy_correction_calculator::material::gas);
    }
    else // xcalo
    {
      // 1 devided by the tracke angle measured relatively to the source foil normal
      double one_over_cos = dist_gas / std::abs(source_vertex_pos_->Y() - calo_vertex_pos_->Y());

      // calculate distance passed in mylar and nylon
      double dist_mylar = one_over_cos * mylar_thickness_;

      // apply energy loss correction (there is no nylon on xwall)
      double E_mylar = Ef_optical + corr_calculator_->Bethe_Bloch_loss(dist_mylar, Ef_optical, energy_correction_calculator::material::mylar);
    	E = E_mylar + corr_calculator_->Bethe_Bloch_loss(dist_gas, E_mylar, energy_correction_calculator::material::gas);
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