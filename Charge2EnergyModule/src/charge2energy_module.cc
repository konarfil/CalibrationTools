#include "charge2energy_module.h"

// Falaise
#include <falaise/snemo/datamodels/precalibrated_data.h>
#include <falaise/snemo/datamodels/calibrated_data.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/datamodels/event.h>
#include <falaise/snemo/datamodels/geomid_utils.h>

// Bayeux
#include <geomtools/geometry_service.h>
#include <datatools/clhep_units.h>

// Other
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

DPP_MODULE_REGISTRATION_IMPLEMENT(charge2energy_module, "Charge2EnergyModule")

charge2energy_module::charge2energy_module()
{

}

charge2energy_module::~charge2energy_module()
{

}

void charge2energy_module::initialize (const datatools::properties & module_properties, datatools::service_manager & services, dpp::module_handle_dict_type &)
{
  event_counter_ = 0;
  geo_manager_ = snemo::service_handle<snemo::geometry_svc>{services};

  double gas_pressure, He_pressure, Et_pressure, Ar_pressure, T_gas;
  std::string calibration_path;

  // read logging priority from the conf file
  if(module_properties.has_key("logging.priority"))
  {
    datatools::logger::priority prio = datatools::logger::get_priority(module_properties.fetch_string("logging.priority"));
    set_logging_priority(prio);
  }
  else
  {
    set_logging_priority(datatools::logger::priority::PRIO_FATAL);
  }

   // read calibration parameter text file path from the conf file
  if(module_properties.has_key("calibration_path"))
    calibration_path = module_properties.fetch_string("calibration_path");
  else
    DT_THROW(std::runtime_error, "Missing calibration parameter file path !");

  // read gas properties from the conf file
  if(module_properties.has_key("gas_pressure"))
    gas_pressure = module_properties.fetch_real("gas_pressure");
  else
    gas_pressure = 0.89;

  if(module_properties.has_key("He_pressure"))
    He_pressure = module_properties.fetch_real("He_pressure");
  else
    He_pressure = 0.955;

  if(module_properties.has_key("Et_pressure"))
    Et_pressure = module_properties.fetch_real("Et_pressure");
  else
    Et_pressure = 0.035;

  if(module_properties.has_key("Ar_pressure"))
    Ar_pressure = module_properties.fetch_real("Ar_pressure");
  else
    Ar_pressure = 0.01;

  if(module_properties.has_key("T_gas"))
    T_gas = module_properties.fetch_real("T_gas");
  else
    T_gas = 298.0;

  if(std::abs(He_pressure + Et_pressure + Ar_pressure - 1.0) > 0.000001)
    DT_LOG_WARNING(get_logging_priority(), "Tracking gas partial pressures do not add up to 1");

  // Initialize the pol3d parameters for MWall 8"
  std::string pol3d_parameters_mwall_8inch_path = "@falaise:snemo/demonstrator/reconstruction/db/fit_parameters_10D_MW_8inch.db";
  if (module_properties.has_key("pol3d_parameters_mwall_8inch_path"))
  {
    pol3d_parameters_mwall_8inch_path = module_properties.fetch_string("pol3d_parameters_mwall_8inch_path");
  }
  datatools::fetch_path_with_env(pol3d_parameters_mwall_8inch_path);
  // Initialize the pol3d parameters for MWall 5"
  std::string pol3d_parameters_mwall_5inch_path = "@falaise:snemo/demonstrator/reconstruction/db/fit_parameters_10D_MW_5inch.db";
  if (module_properties.has_key("pol3d_parameters_mwall_5inch_path"))
  {
    pol3d_parameters_mwall_5inch_path = module_properties.fetch_string("pol3d_parameters_mwall_5inch_path");
  }
  datatools::fetch_path_with_env(pol3d_parameters_mwall_5inch_path);
  // Initialize the pol3d parameters for XWall
  std::string pol3d_parameters_xwall_path = "@falaise:snemo/demonstrator/reconstruction/db/fit_parameters_10D_XW.db";
  if (module_properties.has_key("pol3d_parameters_xwall_path"))
  {
    pol3d_parameters_xwall_path = module_properties.fetch_string("pol3d_parameters_xwall_path");
  }
  datatools::fetch_path_with_env(pol3d_parameters_xwall_path);

  corr_calculator_ = new energy_correction_calculator(gas_pressure, He_pressure, Et_pressure, Ar_pressure,
                                                      T_gas, pol3d_parameters_mwall_8inch_path,
                                                      pol3d_parameters_mwall_5inch_path, pol3d_parameters_xwall_path);

  // initialize calibration parameter vector
  calibration_params_.reserve(number_of_OMs_);
	for (int om = 0; om < number_of_OMs_; om++)
	  calibration_params_.push_back({});

  parse_calibration_params(calibration_path, calibration_params_);

  this->_set_initialized(true);
}

dpp::chain_module::process_status charge2energy_module::process (datatools::things & event)
{
  if(event_counter_ % 10000 == 0)
  {
    DT_LOG_INFORMATION(get_logging_priority(), "Event no. " + std::to_string(event_counter_) + " processed");
  }

  // Skip processing if pCD bank is not present
  if (!event.has("pCD"))
  {
    DT_LOG_WARNING(get_logging_priority(), "======== no pCD bank in event " + std::to_string(event_counter_++) + " ========");
    return dpp::base_module::PROCESS_SUCCESS;
  }

  // Skip processing if PTD bank is not present
  if (!event.has("PTD"))
  {
    DT_LOG_WARNING(get_logging_priority(), "======== no PTD bank in event " + std::to_string(event_counter_++) + " ========");
    return dpp::base_module::PROCESS_SUCCESS;
  }

  // Skip processing if CD bank is not present
  // CD must be present to have track-calo hit association
  if (!event.has("CD"))
  {
    DT_LOG_WARNING(get_logging_priority(), "======== no CD bank in event " + std::to_string(event_counter_++) + " ========");
    return dpp::base_module::PROCESS_SUCCESS;
  }

  const snemo::datamodel::precalibrated_data & pCD = event.get<snemo::datamodel::precalibrated_data>("pCD");
  const snemo::datamodel::particle_track_data & PTD = event.get<snemo::datamodel::particle_track_data>("PTD");

  // calibrate associated calo hits
  for (const datatools::handle<snemo::datamodel::particle_track> & particle_track : PTD.particles())
  {
    for(datatools::handle<snemo::datamodel::calibrated_calorimeter_hit> cd_calo_hit : particle_track->get_associated_calorimeter_hits())
    {
      const datatools::handle<snemo::datamodel::precalibrated_calorimeter_hit> & pcd_calo_hit = pCD.calorimeter_hits()[cd_calo_hit->get_hit_id()];
      calibrate_associated_calo_hit(pcd_calo_hit.get(), particle_track.get(), cd_calo_hit.grab());
    }
  }

  // calibrate isolated calo hits
  for(datatools::handle<snemo::datamodel::calibrated_calorimeter_hit> cd_calo_hit : PTD.isolatedCalorimeters())
  {
    const datatools::handle<snemo::datamodel::precalibrated_calorimeter_hit> & pcd_calo_hit = pCD.calorimeter_hits()[cd_calo_hit->get_hit_id()];
    calibrate_isolated_calo_hit(pcd_calo_hit.get(), cd_calo_hit.grab());
  }
  event_counter_++;

  return dpp::base_module::PROCESS_SUCCESS;
}

void charge2energy_module::parse_calibration_params(std::string database_path, std::vector<std::vector<double>> & params)
{
  std::ifstream database_file (database_path.c_str());
  if(!database_file.is_open())
    DT_THROW(std::runtime_error, "Could not open calibration parameter file !");

  int N_entries = 0;

  std::string line;

  while (std::getline(database_file, line))
  {
    std::stringstream stream (line);
    std::string token;
    std::vector<double> values;

    if (line[0] == '#')
	      continue;

    while (std::getline(stream, token, ';'))
    {
      values.push_back(std::stod(token));
    }

    int om_num = static_cast<int>(values[0]);

    params[om_num] = {values[1], values[2]};
    N_entries++;
  }
  if(N_entries != number_of_OMs_)
    DT_LOG_WARNING(get_logging_priority(), "OMs without calibration parameters will be skipped in energy calculation");
}

// calculate energy for calo hit with associated track
void charge2energy_module::calibrate_associated_calo_hit(const snemo::datamodel::precalibrated_calorimeter_hit & pcd_calo_hit,
                                                         const snemo::datamodel::particle_track & particle_track,
					                                               snemo::datamodel::calibrated_calorimeter_hit & cd_calo_hit)
{
  geomtools::geom_id OM_gid = pcd_calo_hit.get_geom_id();
  int om_num = snemo::datamodel::om_num(OM_gid);
  double charge = -pcd_calo_hit.get_charge() * charge2nVs_;

  // skip if there are no calibration parameters, if the energy is already filled or if the charge is negative
  if(calibration_params_[om_num].size() == 0 || !std::isnan(cd_calo_hit.get_energy()) || charge < 0.0)
    return;

  // find associated calo vertex
  datatools::handle<snemo::datamodel::vertex> calo_vertex;
  for(const datatools::handle<snemo::datamodel::vertex> & vertex : particle_track.get_vertices())
  {
    if(vertex->get_geom_id() == OM_gid)
    {
      calo_vertex = vertex;
    }
  }
  const geomtools::vector_3d & calo_vertex_pos = calo_vertex->get_spot().get_placement().get_translation();

  if(OM_gid.get_type() == 1302) // set the "part" number for main wall modules
  {
    const geomtools::id_mgr & mgr = geo_manager_->get_id_mgr();
    mgr.set(OM_gid, "part", 1);
  }

  // find relative vertex position for optical correction
  const geomtools::mapping   & mapping = geo_manager_->get_mapping();
  const geomtools::geom_info & a_block_ginfo  = mapping.get_geom_info(OM_gid);
  const geomtools::placement & a_block_world_placement = a_block_ginfo.get_world_placement();

  geomtools::placement calo_vertex_pos_OM_placement;
  a_block_world_placement.relocate(calo_vertex_pos, calo_vertex_pos_OM_placement);
  geomtools::vector_3d calo_vertex_pos_OM_vec = calo_vertex_pos_OM_placement.get_translation();

  const double position_x =  (calo_vertex_pos_OM_vec.x());
  const double position_y =  (calo_vertex_pos_OM_vec.y());
  const double position_z = -(calo_vertex_pos_OM_vec.z());

  double position_xyz[3] = {position_y, -position_x, position_z};
  
  // calculate geometrical non-uniformity factor
  double geometrical_factor = corr_calculator_->geometrical_factor(position_xyz[0], position_xyz[1], position_xyz[2], pcd_calo_hit.get_geom_id().get_type());

  // observed energy without corrections
  double Ef = calibration_params_[om_num][0] * charge + calibration_params_[om_num][1];

  // energy with Birks and Cherenkov correction
  double Ef_bc = corr_calculator_->non_linearity_correction(Ef, 1.0);

  // energy with optical correction
  double Ef_optical = corr_calculator_->non_linearity_correction(Ef, geometrical_factor);

  // save different energies into auxiliaries
  cd_calo_hit.grab_auxiliaries().store_with_explicit_unit("Ef", Ef * CLHEP::keV, "Energy without corrections");
  cd_calo_hit.grab_auxiliaries().store_with_explicit_unit("Ef_bc", Ef_bc * CLHEP::keV, "Energy with Birks and Cherenkov corrections");
  cd_calo_hit.grab_auxiliaries().store_with_explicit_unit("Ef_optical", Ef_optical * CLHEP::keV, "Energy optical correction");

  // set energy to the energy without corrections
  cd_calo_hit.set_energy(Ef * CLHEP::keV);

  // find vertex on reference plane
  datatools::handle<snemo::datamodel::vertex> source_vertex;
  for(const datatools::handle<snemo::datamodel::vertex> & vertex : particle_track.get_vertices())
  {
    if(vertex->is_on_reference_source_plane())
    {
      source_vertex = vertex;
    }
  }

  // we only apply energy loss correction if there is a vertex on reference source plane
  // we also skip kinked tracks - to be implemented in the future
  if(source_vertex.has_data() && particle_track.get_trajectory_handle()->get_pattern().number_of_kinks() == 0)
  {
    const geomtools::vector_3d & source_vertex_pos = source_vertex->get_spot().get_placement().get_translation();
    double dist_gas = std::sqrt((calo_vertex_pos.getX() - source_vertex_pos.getX())*(calo_vertex_pos.getX() - source_vertex_pos.getX())
                              + (calo_vertex_pos.getY() - source_vertex_pos.getY())*(calo_vertex_pos.getY() - source_vertex_pos.getY())
                              + (calo_vertex_pos.getZ() - source_vertex_pos.getZ())*(calo_vertex_pos.getZ() - source_vertex_pos.getZ()));

    // 1 divided by the track angle measured relatively to the source foil normal
    double one_over_cos = dist_gas / std::abs(source_vertex_pos.getX() - calo_vertex_pos.getX());

    // calculate distance passed in mylar and nylon
    double dist_mylar = one_over_cos * mylar_thickness_;
    double dist_nylon = one_over_cos * nylon_thickness_;

    
    double Ef_mylar = Ef_optical + corr_calculator_->Bethe_Bloch_loss(dist_mylar, Ef_optical, energy_correction_calculator::material::mylar);
    double Ef_nylon = Ef_mylar + corr_calculator_->Bethe_Bloch_loss(dist_nylon, Ef_mylar, energy_correction_calculator::material::nylon);

    // energy with optical correction + energy loss correction
    double Ef_gas   = Ef_nylon + corr_calculator_->Bethe_Bloch_loss(dist_gas, Ef_nylon, energy_correction_calculator::material::gas);
    cd_calo_hit.grab_auxiliaries().store_with_explicit_unit("Ef_optical_loss", Ef_gas * CLHEP::keV, "Energy with both optical and energy loss corrections");
  }
}

// calculate energy for calo hit without associated track
void charge2energy_module::calibrate_isolated_calo_hit(const snemo::datamodel::precalibrated_calorimeter_hit & pcd_calo_hit,
					                                             snemo::datamodel::calibrated_calorimeter_hit & cd_calo_hit)
{
  int om_num = snemo::datamodel::om_num(pcd_calo_hit.get_geom_id());
  double charge = -pcd_calo_hit.get_charge() * charge2nVs_;

  // skip if there are no calibration parameters, if the energy is already filled or if the charge is negative
  if(calibration_params_[om_num].size() == 0.0 || !std::isnan(cd_calo_hit.get_energy()) || charge < 0.0)
    return;

  // energy without energy corrections
  double Ef = calibration_params_[om_num][0] * charge + calibration_params_[om_num][1];

  // energy with Birks and Cherenkov correction
  double Ef_bc = corr_calculator_->non_linearity_correction(Ef, 1.0);

  cd_calo_hit.set_energy(Ef * CLHEP::keV);

  cd_calo_hit.grab_auxiliaries().store_with_explicit_unit("Ef", Ef * CLHEP::keV, "Energy without corrections");
  cd_calo_hit.grab_auxiliaries().store_with_explicit_unit("Ef_bc", Ef_bc * CLHEP::keV, "Energy with Birks and Cherenkov corrections");
}