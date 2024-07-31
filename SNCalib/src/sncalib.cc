#include "calibration_parameter_finder.h"
#include "energy_correction_calculator.h"
#include "calib_info.h"

// ROOT
#include <TFile.h>
#include <TTree.h>

// Other
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

// load parameters from the configuration file
std::map<std::string, std::string> read_config_file(const std::string& filename)
{
  std::map<std::string, std::string> config;
  std::ifstream file(filename);
  std::string line;
  while (std::getline(file, line))
  {
    std::istringstream is_line(line);
    std::string key;
    if (std::getline(is_line, key, '=')) 
    {
      std::string value;
      if (std::getline(is_line, value)) 
      {
        config[key] = value;
      }
    }
  }
  return config;
}

int main(int argc, char *argv[])
{
  gErrorIgnoreLevel = kError; // lower ROOT verbosity

  std::string input_file_path = "";
  std::string output_file_path = "";
  std::string config_file_path = "";
  bool save_fitted_spectra = false;
  bool verbose = false;

  for (int iarg=1; iarg<argc; ++iarg)
  {
    std::string arg (argv[iarg]);
    if (arg[0] == '-')
    {
      if (arg=="-i" || arg=="--input")
      {
        input_file_path = std::string(argv[++iarg]);
      }
      else if (arg=="-h" || arg=="--help")
      {
        std::cout << "-h [--help]\tprint help messenge" << std::endl;
        std::cout << "-i [--input]\tinput root file with calibration data for individual OMs" << std::endl;
        std::cout << "-o [--output]\toutput csv file to save calibration parameters into" << std::endl;
        std::cout << "-p [--config]\tconfiguration file with parameters of the calibration algorithm" << std::endl;
        std::cout << "-s [--spectra]\tflag to save fitted spectra as png files" << std::endl;
        std::cout << "-V [--verbose]\tflag to print number of processed OM" << std::endl;
        return 0;
      }
      else if (arg=="-o" || arg=="--output")
      {
        output_file_path = std::string(argv[++iarg]);
      }
      else if (arg=="-p" || arg=="--config")
      {
        config_file_path = std::string(argv[++iarg]);
      }
      else if (arg=="-s" || arg=="--spectra")
      {
        save_fitted_spectra = true;
      }
      else if (arg=="-V" || arg=="--verbose")
      {
        verbose = true;
      }
      else
      {
        std::cerr << "unkown option " << arg << std::endl;
      }
    }
  }

  if (input_file_path.empty())
  {
    std::cerr << "missing input file !" << std::endl;
    return 1;
  }

  if (output_file_path.empty())
  {
    std::cerr << "missing output file !" << std::endl;
    return 1;
  }

  if (config_file_path.empty())
  {
    std::cerr << "missing configuration file !" << std::endl;
    return 1;
  }
  std::map<std::string, std::string> config = read_config_file(config_file_path);
  int min_hits = std::stoi(config["min_hits"]);
  double minimization_threshold = std::stod(config["minimization_threshold"]);

  // geometrical non-uniformity correction paths
  std::string falaise_resources_path = Falaise_RESOURCE_DIR;
  std::string pol3d_parameters_mwall_8inch_path = falaise_resources_path + "/snemo/demonstrator/reconstruction/db/fit_parameters_10D_MW_8inch.db";
  std::string pol3d_parameters_mwall_5inch_path = falaise_resources_path + "/snemo/demonstrator/reconstruction/db/fit_parameters_10D_MW_5inch.db";
  std::string pol3d_parameters_xwall_path = falaise_resources_path + "/snemo/demonstrator/reconstruction/db/fit_parameters_10D_XW.db";

  energy_correction_calculator* corr_calculator = new energy_correction_calculator(std::stod(config["gas_pressure"]), std::stod(config["He_pressure"]),
                                                      std::stod(config["Et_pressure"]), std::stod(config["Ar_pressure"]),
                                                      std::stod(config["T_gas"]), pol3d_parameters_mwall_8inch_path,
                                                      pol3d_parameters_mwall_5inch_path, pol3d_parameters_xwall_path);

  TFile* calib_data_file = new TFile(input_file_path.c_str());
  std::ofstream calib_param_file;
  calib_param_file.open (output_file_path);
  calib_param_file << "#OM_number;a;b;chi2_A;chi2_B;loss" << std::endl;
  for(auto keyObj : *calib_data_file->GetListOfKeys())
  {
    TKey* key = (TKey*)keyObj;
    std::string om_title = key->GetName();
    TTree* OM_data = (TTree*)calib_data_file->Get(om_title.c_str());
    if(OM_data->GetEntries() < min_hits)
    {
      continue;
    }

    calibration_parameter_finder finder = calibration_parameter_finder(OM_data, std::stod(config["minimization_threshold"]), 
                                                                       std::stoi(config["max_iterations"]), save_fitted_spectra, corr_calculator);
    calib_info best_calib = finder.find_calibration_parameters();

    calib_param_file << om_title << ";" << best_calib.a << ";" << best_calib.b 
                     << ";" << best_calib.chi2_NDF_A << ";" << best_calib.chi2_NDF_B
                     << ";" << best_calib.loss << std::endl;

    if(verbose) std::cout << "OM " << om_title << " processed" << std::endl;
    if(best_calib.loss > minimization_threshold)
      std::cout << "Warning: OM " << om_title << " did not converge. Lowest loss function value found (" << best_calib.loss << ") is higher than chosen threshold (" << minimization_threshold << ")" << std::endl;
  }
  calib_param_file.close();
  delete calib_data_file;
  delete corr_calculator;

  return 0;
}

