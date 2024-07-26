#ifndef CALIBRATION_PARAMETER_FINDER_H
#define CALIBRATION_PARAMETER_FINDER_H

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
    calibration_parameter_finder(TTree* OM_data, double gas_pressure, double He_pressure, double Et_pressure, 
                                double Ar_pressure, double T_gas, double minimization_threshold, int max_iterations,
                                bool save_fitted_spectra);

    // Destructor
    ~calibration_parameter_finder();

    // calculate calibration parameters of the OM and assign them to a and b
    calib_info find_calibration_parameters();

  private:
    // fitting function for the first Bi207 peak (sum of three gaussians)
    static double tripleGaussA(double* x,  double* par);

    // fitting function for the second Bi207 peak (sum of three gaussians)
    static double tripleGaussB(double* x,  double* par);

    // parse parameters for the non-uniformity correction
    std::vector<double> parse_pol3d_parameters(const std::string & parameters_path);

    // calculate energy losses in material using the Bethe-Bloch formula for electrons
    double Bethe_Bloch_loss(double d, double E, double Z_A, double I, double RHO);

    // based on observed energy and geometrical correction factor 
    // calculates energy corrected by the non-linearity effects
    double non_linearity_correction(double Ef, double geometrical_factor);

    // a function of calibration parameters to be minimized
    double loss_function(double a, double b);

    // intensities of Bi207 internal conversion electrons
    static constexpr double I_482keV_  = 1.52;
    static constexpr double I_555keV_  = 0.44;
    static constexpr double I_567keV_  = 0.15;
    static constexpr double I_976keV_  = 7.03;
    static constexpr double I_1049keV_ = 1.84;
    static constexpr double I_1061keV_ = 0.54;

    static constexpr double E_0_ = 511.0; // rest energy of electron

    // thicknesses of material layers in mm
    static constexpr double mylar_thickness_ = 0.012;
    static constexpr double nylon_thickness_ = 0.025;
    static constexpr double gas_thickness_   = 435.0;

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

    // coeffitiens of the minimization algorithm
    static constexpr double ALPHA_ = 1.0;   // Reflection coefficient
    static constexpr double GAMMA_ = 2.0;   // Expansion coefficient
    static constexpr double RHO_ = 0.5;     // Contraction coefficient
    static constexpr double SIGMA_ = 0.5;   // Shrink coefficient

    // gas properites for energy loss calculation using the Bethe-Bloch formula
    double Z_A_gas_;
    double I_gas_;
    double rho_gas_;

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