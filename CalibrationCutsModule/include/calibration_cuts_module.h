#ifndef CALIBRATION_CUTS_MODULE_H
#define CALIBRATION_CUTS_MODULE_H

// Falaise
#include <falaise/snemo/services/service_handle.h>
#include <falaise/snemo/services/geometry.h>
#include <falaise/snemo/datamodels/particle_track.h>

// Bayeux
#include <bayeux/dpp/chain_module.h>
#include <geomtools/manager.h>

// ROOT
#include <TTree.h>
#include <TVector3.h>
#include <TFile.h>
#include <TH2F.h>

// Other
#include <string>
#include <map>
#include <memory>

class calibration_cuts_module : public dpp::chain_module
{
  public:
    // Constructor
    calibration_cuts_module();
  
    // Destructor
    virtual ~calibration_cuts_module();
  
    // Initialisation function
    virtual void initialize (const datatools::properties &,
                             datatools::service_manager &,
  			                     dpp::module_handle_dict_type &);
  
    // Event processing function
    dpp::chain_module::process_status process (datatools::things & event);
    
  private:

    // check if the particle track has exactly one associated OM hit - if so, return the vertex associated with the hit
    datatools::handle<snemo::datamodel::vertex> track_has_one_assoc_calo(const datatools::handle<snemo::datamodel::particle_track> & particle);

    // check if the particle track has a vertex near a calibration source - if so, return the vertex
    datatools::handle<snemo::datamodel::vertex> vertex_close_to_a_calib_source(const datatools::handle<snemo::datamodel::particle_track> & particle);

    // based on OM number return a string with GID
    std::string om_num_to_gid(int om_num);

    std::string output_path_; // path of a file to save extracted information into

    int event_counter_;
    snemo::service_handle<snemo::geometry_svc> geo_manager_{};

    static const int number_of_OMs_ = 648; // total number of OMs in SuperNEMO without gveto
    TTree* OM_data_[number_of_OMs_]; // an array of TTrees to save calibration information into

    static constexpr double charge2nVs_ = 1e6; // multiplicative constant to transform charge to units nV*s

    static constexpr double charge_threshold_ = 0.5; // minumum charge to be saved

    // number of calibration sources
    static const int calib_source_rows_ = 7;
    static const int calib_source_columns_ = 6;

    // dimensions of the ellipse for calibration source vertex cut
    double source_cut_ellipse_Y_;
    double source_cut_ellipse_Z_;

    // Y and Z coorinates of centres of calibration sources
    double calib_source_Y_[calib_source_rows_][calib_source_columns_];
    double calib_source_Z_[calib_source_rows_][calib_source_columns_];

    // variables for saving data into TTree branches
    float charge_; // charge measured by the OM
    TVector3 source_vertex_pos_; // position of the calibration source vertex
    TVector3 calo_vertex_pos_; // position of the calorimeter vertex
    TVector3 calo_vertex_pos_OM_; // position of the calorimeter vertex in coordinate 
                                 // system of the OM (needed for the optical correction)

    // file to save extracted information into
    TFile* save_file_;

    int far_from_source_ = 0;
  
    // Macro to register the module
    DPP_MODULE_REGISTRATION_INTERFACE(calibration_cuts_module);
};
#endif