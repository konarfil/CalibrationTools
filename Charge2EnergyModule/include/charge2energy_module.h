#ifndef CHARGE2ENERGY_MODULE_H
#define CHARGE2ENERGY_MODULE_H

#include "energy_correction_calculator.h"

// Falaise
#include <falaise/snemo/datamodels/precalibrated_calorimeter_hit.h>
#include <falaise/snemo/datamodels/calibrated_calorimeter_hit.h>
#include <falaise/snemo/datamodels/particle_track.h>
#include <falaise/snemo/services/service_handle.h>
#include <falaise/snemo/services/geometry.h>

// Bayeux
#include <bayeux/dpp/chain_module.h>
#include <geomtools/manager.h>

// Other
#include <vector>

class charge2energy_module : public dpp::chain_module
{
  public:
    // Constructor
    charge2energy_module();
  
    // Destructor
    virtual ~charge2energy_module();
  
    // Initialisation function
    virtual void initialize (const datatools::properties &,
                             datatools::service_manager &,
  			                     dpp::module_handle_dict_type &);
  
    // Event processing function
    dpp::chain_module::process_status process (datatools::things & event);
    
  private:
    void parse_calibration_params(std::string database_path_, std::vector<std::vector<double>> & params_);

    // calculate energy for calo hit with associated track
    void calibrate_associated_calo_hit(const snemo::datamodel::precalibrated_calorimeter_hit & pcd_calo_hit,
                                       const snemo::datamodel::particle_track & particle_track,
					                             snemo::datamodel::calibrated_calorimeter_hit & cd_calo_hit);

    // calculate energy for calo hit without associated track
    void calibrate_isolated_calo_hit(const snemo::datamodel::precalibrated_calorimeter_hit & pcd_calo_hit,
					                           snemo::datamodel::calibrated_calorimeter_hit & cd_calo_hit);

    energy_correction_calculator* corr_calculator_;

    snemo::service_handle<snemo::geometry_svc> geo_manager_{};

    static const int number_of_OMs_ = 712; // total number of OMs in SuperNEMO
    static constexpr double charge2nVs_ = 1e6; // multiplicative constant to transform charge to units nV*s

    // thicknesses of material layers in mm
    static constexpr double mylar_thickness_ = 0.012;
    static constexpr double nylon_thickness_ = 0.025;

    int event_counter_;

    // calibration parameters of individual OMs
    std::vector<std::vector<double>> calibration_params_;
  
    // Macro to register the module
    DPP_MODULE_REGISTRATION_INTERFACE(charge2energy_module);
};
#endif