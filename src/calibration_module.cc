#include "calibration_module.h"

// Falaise
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/datamodels/precalibrated_data.h>
#include <falaise/snemo/datamodels/geomid_utils.h>

// Bayeux
#include <geomtools/geometry_service.h>
#include <datatools/clhep_units.h>

// Macro to add the module in the global register of data processing modules:
// The module defined by this class 'falaise_skeleton_module_ptd' will be registered
// with the label ID 'FalaiseSkeletonModule_PTD' (to use in pipeline configuration file)
DPP_MODULE_REGISTRATION_IMPLEMENT(calibration_module, "CalibrationModule")

calibration_module::calibration_module()
{
  std::cout << "calibration_module::calibration_module() called" << std::endl;
}


calibration_module::~calibration_module()
{
  for(auto item : (*OM_data_))
  {
    for(int i = 0;i < item.second.charge.size();i++)
    {
      std::cout << item.second.charge[i] << " " 
                << item.second.vertex_on_OM_X[i] << " " 
                << item.second.vertex_on_OM_Y[i] << " " 
                << item.second.one_over_cos[i] << std::endl;
    }
  }

  std::cout << "calibration_module::~calibration_module() called" << std::endl;
}

void calibration_module::initialize (const datatools::properties & module_properties, datatools::service_manager & services, dpp::module_handle_dict_type &)
{
  std::cout << "calibration_module::initialize() called" << std::endl;

  event_counter_ = 0;

  geo_manager_ = snemo::service_handle<snemo::geometry_svc>{services};
  OM_data_ = std::make_shared<std::map<int,OM_data>>();

  // read calibration source cut parameters from conf file
  if ( module_properties.has_key("source_cut_ellipse_Y"))
    source_cut_ellipse_Y_ = module_properties.fetch_real("source_cut_ellipse_Y") / CLHEP::mm;
  else
    source_cut_ellipse_Y_ = 25.0 * CLHEP::mm;

  if ( module_properties.has_key("source_cut_ellipse_Z"))
    source_cut_ellipse_Z_ = module_properties.fetch_real("source_cut_ellipse_Z") / CLHEP::mm;
  else
    source_cut_ellipse_Z_ = 30.0 * CLHEP::mm;


  calib_cuts_driver_ = new calibration_cuts_driver(OM_data_, geo_manager_, source_cut_ellipse_Y_, source_cut_ellipse_Z_);

  this->_set_initialized(true);
}

dpp::chain_module::process_status calibration_module::process (datatools::things & event)
{
  // Skip processing if PTD bank is not present
  if (!event.has("PTD"))
  {
    std::cout << "======== no PTD bank in event " << event_counter_++ << " ========" << std::endl;
    return dpp::base_module::PROCESS_SUCCESS;
  }

  // Skip processing if pCD bank is not present
  if (!event.has("pCD"))
  {
    std::cout << "======== no pCD bank in event " << event_counter_++ << " ========" << std::endl;
    return dpp::base_module::PROCESS_SUCCESS;
  }

  // Retrieve data banks
  const snemo::datamodel::particle_track_data & PTD = event.get<snemo::datamodel::particle_track_data>("PTD");
  const snemo::datamodel::precalibrated_data & pCD = event.get<snemo::datamodel::precalibrated_data>("pCD");

  // apply data cuts and save information into OM_data_
  calib_cuts_driver_->process_event(PTD, pCD);

  event_counter_++;

  return dpp::base_module::PROCESS_SUCCESS;
}