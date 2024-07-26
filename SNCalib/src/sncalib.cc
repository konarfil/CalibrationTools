#include "calibration_parameter_finder.h"
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

  TFile* calib_data_file = new TFile(input_file_path.c_str());
  std::ofstream calib_param_file;
  calib_param_file.open (output_file_path);
  for(auto keyObj : *calib_data_file->GetListOfKeys())
  {
    TKey* key = (TKey*)keyObj;
    std::string om_title = key->GetName();
    TTree* OM_data = (TTree*)calib_data_file->Get(om_title.c_str());
    if(OM_data->GetEntries() < min_hits)
    {
      continue;
    }

    calibration_parameter_finder finder = calibration_parameter_finder(OM_data, std::stod(config["gas_pressure"]), std::stod(config["He_pressure"]),
                                                                       std::stod(config["Et_pressure"]), std::stod(config["Ar_pressure"]), std::stod(config["T_gas"]),
                                                                       std::stod(config["minimization_threshold"]), std::stoi(config["max_iterations"]),
                                                                       save_fitted_spectra);
    calib_info best_calib = finder.find_calibration_parameters();

    calib_param_file << om_title << ";" << best_calib.a << ";" << best_calib.b 
                     << ";" << best_calib.chi2_NDF_A << ";" << best_calib.chi2_NDF_B
                     << ";" << best_calib.loss << std::endl;
    std::cout << om_title << " processed" << std::endl;
  }
  calib_param_file.close();
  delete calib_data_file;

  return 0;
}

