#ifndef CALIBRATION_PARAMETER_FINDER_H
#define CALIBRATION_PARAMETER_FINDER_H

#include "energy_correction_calculator.h"

// Falaise
#include <falaise/snemo/services/service_handle.h>
#include <falaise/snemo/services/geometry.h>

// Bayeux
#include <geomtools/manager.h>
#include <geomtools/geometry_service.h>

// ROOT
#include <TTree.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TF1.h>

// Other
#include <vector>
#include <string>

#include "calib_info.h"

class calibration_parameter_finder
{
  public:
    // Constructor
    calibration_parameter_finder(TTree* OM_data, double minimization_threshold, int max_iterations,
                                bool save_fitted_spectra, energy_correction_calculator* corr_calculator);

    // Destructor
    ~calibration_parameter_finder();

    // calculate calibration parameters of the OM and assign them to a and b
    calib_info find_calibration_parameters();

  private:
    // fitting function for the first Bi207 peak (sum of three gaussians)
    static double tripleGaussA(double* x,  double* par);

    // fitting function for the second Bi207 peak (sum of three gaussians)
    static double tripleGaussB(double* x,  double* par);

    // a function of calibration parameters to be minimized
    double loss_function(double a, double b);

    // coeffitiens of the minimization algorithm
    static constexpr double ALPHA_ = 1.0;   // Reflection coefficient
    static constexpr double GAMMA_ = 2.0;   // Expansion coefficient
    static constexpr double RHO_ = 0.5;     // Contraction coefficient
    static constexpr double SIGMA_ = 0.5;   // Shrink coefficient

    // thicknesses of material layers in mm
    static constexpr double mylar_thickness_ = 0.012;
    static constexpr double nylon_thickness_ = 0.025;

    // intensities of Bi207 internal conversion electrons
    static constexpr double I_482keV_  = 1.52;
    static constexpr double I_555keV_  = 0.44;
    static constexpr double I_567keV_  = 0.15;
    static constexpr double I_976keV_  = 7.03;
    static constexpr double I_1049keV_ = 1.84;
    static constexpr double I_1061keV_ = 0.54;

    energy_correction_calculator* corr_calculator_;

    TTree* OM_data_; // data corresponding to indivitual OMs

    // variables for OM_data_ branches reading
    float charge_;
    TVector3* source_vertex_pos_;
    TVector3* calo_vertex_pos_;
    TVector3* calo_vertex_pos_OM_;

    double minimization_threshold_; // threshold of the loss function to end the minimization
    int max_iterations_; // maximum number of iterations for the minimization algorithm
    bool save_fitted_spectra_; // save fitted calibration spectra as pngs ?
    std::string fitted_spectra_dir_; // directory to save pngs of fitted spectra into

    calib_info best_calib_; // information about the best found calibration

    TH1F min_loss_histo_; // calibration spectrum corresponding to the lowest found loss
    TF1  min_loss_gaussA_; // fitting function of the 482 keV peak corresponding to the lowest found loss
    TF1  min_loss_gaussB_; // fitting function of the 976 keV peak corresponding to the lowest found loss

    // initial points of the initial algorithm
    static const int N_minimization_points_ = 3;
    static const int N_calibration_parameters_ = 2;
      
    std::vector<double> uniformity_correction_parameters_mwall_8inch_{1,1}; // Polynomial parameters for the uniformity correction for MWall 8"
    std::vector<double> uniformity_correction_parameters_mwall_5inch_{1,1}; // Polynomial parameters for the uniformity correction for MWall 5"
    std::vector<double> uniformity_correction_parameters_xwall_{1,1};       // Polynomial parameters for the uniformity correction for XWall
};

#endif