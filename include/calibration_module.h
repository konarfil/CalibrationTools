#ifndef CALIBRATION_MODULE_H
#define CALIBRATION_MODULE_H

#include "OM_data.h"
#include "calibration_cuts_driver.h"

// Falaise
#include <falaise/snemo/services/service_handle.h>
#include <falaise/snemo/services/geometry.h>

// Bayeux
#include <bayeux/dpp/chain_module.h>
#include <geomtools/manager.h>

// Other
#include <string>
#include <map>
#include <memory>
class calibration_module : public dpp::chain_module
{
  public:
    // Constructor
    calibration_module();
  
    // Destructor
    virtual ~calibration_module();
  
    // Initialisation function
    virtual void initialize (const datatools::properties &,
                             datatools::service_manager &,
  			                     dpp::module_handle_dict_type &);
  
    // Event processing function
    dpp::chain_module::process_status process (datatools::things & event);
    
  private:
    calibration_cuts_driver* calib_cuts_driver_;

    // dimensions of the ellipse for calibration source vertex cut
    double source_cut_ellipse_Y_;
    double source_cut_ellipse_Z_;

    int event_counter_;
    std::shared_ptr<std::map<int, OM_data>> OM_data_; // information about individual OMs, key = OM number

    snemo::service_handle<snemo::geometry_svc> geo_manager_{};
  
    // Macro to register the module
    DPP_MODULE_REGISTRATION_INTERFACE(calibration_module);
};
#endif