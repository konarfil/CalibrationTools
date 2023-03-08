const double I_976keV  = 7.03;
const double I_1049keV = 1.84;
const double I_1061keV = 0.54;

double tripleGaussB(double* x,  double* par);

int fit_spectra(string input_folder, string output_folder)
{	
	gROOT->SetBatch(kTRUE);
	gErrorIgnoreLevel = 6001; // to ignore warnings about empty histograms
	
	// load directory with OM spectra
	TSystemDirectory dir(input_folder.c_str(), input_folder.c_str());
	TList *files = dir.GetListOfFiles();
	if (files)
	{
		TSystemFile *file;
		string fname;
		TIter next(files);
		
		// skip "." and ".."
		next();
		next();
		
		
		ofstream param_file((output_folder + "calibration_params.txt"));
		while ((file=(TSystemFile*)next()))
		{
			fname = file->GetName();
			TFile* file = new TFile((input_folder + fname).c_str());
			
			// look for spectra with second peak at 976 keV
			double min_k = 0.0;
			double min_error = 1e10;
			double min_sigma = 0.0;
			TH1F   min_histo;
			TF1    min_fit;
			
			for(auto keyObj : *file->GetListOfKeys())
			{
				TKey* key = (TKey*)keyObj;
				TH1F* histo = (TH1F*)file->Get(key->GetName());
				
				// get histogram atributes
				double histo_max   = histo->GetMaximum();
				double bin_width   = histo->GetBinWidth(1);
				double X_min       = histo->GetXaxis()->GetXmin();
				int    N_bins      = histo->GetNbinsX();
				
				// determine fitting range based on running average of the spectrum
				TH1F* running_average = new TH1F("avg", "avg", N_bins, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
				for(int i = 3;i < N_bins - 3;i++)
				{
					double sum = 0.0;
					for(int j = 0;j < 5;j++)
					{
						sum += histo->GetBinContent(i + j - 2);
					}
					running_average->SetBinContent(i, sum / 5.0);
				}
				
				double max_bin = running_average->GetMaximumBin();
				double max     = running_average->GetBinContent(max_bin);
				
				double fit_min, fit_max;
				
				for(int i = max_bin;i > 1;i--)
				{
					if(running_average->GetBinContent(i) < 0.3 * max)
					{
						fit_min = X_min + i * bin_width;
						break;
					}
				}
				
				for(int i = max_bin;i < N_bins;i++)
				{
					if(running_average->GetBinContent(i) < 0.15 * max)
					{
						fit_max = X_min + i * bin_width;
						break;
					}
				}
				
				// initialize fitting function
				TF1* tripGFitB = new TF1("TripleGB", tripleGaussB, fit_min, fit_max, 3);
				tripGFitB->SetNpx(10000);
				tripGFitB->SetParNames("N1_B",  "mu1_B", "sig1_B");
				
				// set initial fitting parameters (rough estimate)
				tripGFitB->SetParameters(histo_max * 360.0, (fit_min + fit_max) / 2.0, 0.8 * (fit_max - fit_min) / 2.0);
				
				histo->Fit(tripGFitB, "QR");
				
				if(abs((tripGFitB->GetParameter(1) - 976.0)) < min_error)
				{
					min_error = abs(tripGFitB->GetParameter(1) - 976.0);
					string histo_name = key->GetName();
					min_k     = stod(histo_name.substr(histo_name.find("=") + 1, histo_name.length()));
					min_histo = *histo;
					min_fit   = *tripGFitB;
					min_sigma = tripGFitB->GetParameter(2);
				}
				delete histo;
				delete tripGFitB;
				delete running_average;
			}
			// save the best fit and the best calibration parameter
			cout << fname.substr(0, fname.length() - 5) << " finished !" << endl;
			param_file << fname.substr(0, fname.length() - 5) << ";" << min_k << endl;
			TCanvas* c = new TCanvas("c", "c", 1280, 720);
			string histo_title = min_histo.GetTitle();
			min_histo.SetTitle((histo_title + "_#sigma=" + to_string(min_sigma)).c_str());
			min_histo.Draw();
			min_fit.Draw("same");
			c->SaveAs((output_folder + fname.substr(0, fname.length() - 5) + ".png").c_str());
			delete c;
		}
		param_file.close();
	}
	else
	{
		cout << "*** invalid input folder !" << endl;
		return 1;
	}
	return 0;
}

// sum of three gaussian functions used to fit the 976 keV peak
double tripleGaussB(double* x,  double* par)
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
	double n2   = (I_1049keV / I_976keV) * n1;
	double mu2  = (1049.0 / 976.0) * mu1;
	double sig2 = sqrt(1049.0 / 976.0) * sig1;
	double n3   = (I_1061keV / I_976keV) * n1;
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
