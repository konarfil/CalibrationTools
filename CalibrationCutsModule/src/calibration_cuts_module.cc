#include "calibration_cuts_module.h"

// Falaise
#include <falaise/snemo/services/service_handle.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/datamodels/precalibrated_data.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/snemo/datamodels/geomid_utils.h>
#include <falaise/snemo/datamodels/vertex_utils.h>

// Bayeux
#include <geomtools/geometry_service.h>
#include <datatools/clhep_units.h>

DPP_MODULE_REGISTRATION_IMPLEMENT(calibration_cuts_module, "CalibrationCutsModule")

calibration_cuts_module::calibration_cuts_module()
{
  
}

calibration_cuts_module::~calibration_cuts_module()
{
  save_file_->cd();
  for(int i = 0;i < number_of_OMs_no_gveto_;i++)
  {
    OM_data_[i]->Write();
    delete OM_data_[i];
  }
  delete save_file_;
}

void calibration_cuts_module::initialize (const datatools::properties & module_properties, datatools::service_manager & services, dpp::module_handle_dict_type & module_dict)
{
  event_counter_ = 0;

  geo_manager_ = snemo::service_handle<snemo::geometry_svc>{services};

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

  // read calibration source cut parameters from conf file
  if(module_properties.has_key("source_cut_ellipse_Y"))
    source_cut_ellipse_Y_ = module_properties.fetch_real("source_cut_ellipse_Y") / CLHEP::mm;
  else
    source_cut_ellipse_Y_ = 25.0 * CLHEP::mm;

  if(module_properties.has_key("source_cut_ellipse_Z"))
    source_cut_ellipse_Z_ = module_properties.fetch_real("source_cut_ellipse_Z") / CLHEP::mm;
  else
    source_cut_ellipse_Z_ = 30.0 * CLHEP::mm;

  // read output file path from the conf file
  if(module_properties.has_key("output_path"))
    output_path_ = module_properties.fetch_string("output_path");
  else
    output_path_ = "extracted_data.root";

  // read calibration source position file path from the conf file
  if(module_properties.has_key("source_pos_path"))
    source_pos_path_ = module_properties.fetch_string("source_pos_path");
  else
    source_pos_path_ = "";

  if(source_pos_path_ == "")
  {
    // iterate through all calibration sources and save their Y and Z positions into arrays
    for(int i = 0;i < calib_source_rows_;i++)
    {
      for(int j = 0;j < calib_source_columns_;j++)
      {
        const geomtools::id_mgr & mgr = geo_manager_->get_id_mgr();
        geomtools::geom_id calib_spot_id;
        mgr.make_id("source_calibration_spot", calib_spot_id);
        mgr.set(calib_spot_id, "module", 0);
        mgr.set(calib_spot_id, "track", j);
        mgr.set(calib_spot_id, "position", i);
      
        const geomtools::mapping & mapping = geo_manager_->get_mapping();
        const geomtools::geom_info & calib_spot_ginfo = mapping.get_geom_info(calib_spot_id);
        const geomtools::placement & calib_spot_placement = calib_spot_ginfo.get_world_placement();
        const geomtools::vector_3d & calib_spot_pos = calib_spot_placement.get_translation();

        calib_source_Y_[i][j] = calib_spot_pos.getY();
        calib_source_Z_[i][j] = calib_spot_pos.getZ();
      }
    }
  }
  else
  {	
  	std::ifstream source_positions_file(source_pos_path_);
    DT_THROW_IF(source_positions_file.fail(), std::logic_error, source_pos_path_ + " does not exist !");

    int i = 0;
    std::string line;

  	while (getline(source_positions_file, line)) 
  	{
  		calib_source_Y_[i / 6][i % 6] = std::stod(line.substr(0, line.find(";") + 1));
  		calib_source_Z_[i / 6][i % 6] = std::stod(line.substr(line.find(";") + 1, line.length()));

  		i++;
  	}
  	source_positions_file.close(); 
  }

  save_file_ = new TFile(output_path_.c_str(), "RECREATE");

  source_vertex_pos_ = TVector3();
  calo_vertex_pos_ = TVector3();
  calo_vertex_pos_OM_ = TVector3();

  // initialize ttrees to save extracted OM data into
  for(int i = 0;i < number_of_OMs_no_gveto_;i++)
  {
    std::string om_num = std::to_string(i);

    OM_data_[i] = new TTree(om_num.c_str(), om_num.c_str());
    OM_data_[i]->Branch("charge", &charge_);
    OM_data_[i]->Branch("source_vertex_pos", &source_vertex_pos_);
    OM_data_[i]->Branch("calo_vertex_pos", &calo_vertex_pos_);
    OM_data_[i]->Branch("calo_vertex_pos_OM", &calo_vertex_pos_OM_);
  }

  this->_set_initialized(true);
}

dpp::chain_module::process_status calibration_cuts_module::process (datatools::things & event)
{
  if(event_counter_ % 10000 == 0)
  {
    DT_LOG_INFORMATION(get_logging_priority(), "Event no. " + std::to_string(event_counter_) + " processed");
  }

  // Skip processing if PTD bank is not present
  if (!event.has("PTD"))
  {
    DT_LOG_WARNING(get_logging_priority(), "======== no PTD bank in event " + std::to_string(event_counter_++) + " ========");
    return dpp::base_module::PROCESS_SUCCESS;
  }

  // Skip processing if pCD bank is not present
  if (!event.has("pCD"))
  {
    DT_LOG_WARNING(get_logging_priority(), "======== no pCD bank in event " + std::to_string(event_counter_++) + " ========");
    return dpp::base_module::PROCESS_SUCCESS;
  }

  // Retrieve data banks
  const snemo::datamodel::particle_track_data & PTD = event.get<snemo::datamodel::particle_track_data>("PTD");
  const snemo::datamodel::precalibrated_data & pCD = event.get<snemo::datamodel::precalibrated_data>("pCD");

  // iterate through all particles and find those coming from a calibration source and hitting an OM
  for (const datatools::handle<snemo::datamodel::particle_track> & particle : PTD.particles())
	{
    datatools::handle<snemo::datamodel::vertex> source_vertex = vertex_close_to_a_calib_source(particle);
    if(!source_vertex.has_data()) // there is no vertex close to a calibration source
    {
      continue;
    }

    datatools::handle<snemo::datamodel::vertex> calo_vertex = track_has_one_assoc_calo(particle);
    if(!calo_vertex.has_data()) // track does not have exactly one associated OM hit
    {
      continue;
    }

    // we only take tracks without kinks
    // in the future kinks could be included but it would require modification of the energy loss calculation
    // (add all track points to OM_data and calculate losses for individual segments)
    if(particle->get_trajectory_handle()->get_pattern().number_of_kinks() != 0)
    {
      continue;
    }

    // find the associated OM in pCD and save measured charge and its gid
    float charge = -1.0;
    geomtools::geom_id assoc_gid;
    for(const datatools::handle<snemo::datamodel::precalibrated_calorimeter_hit> & hit : pCD.calorimeter_hits())
    {
      if(hit->get_geom_id() == calo_vertex->get_geom_id())
      {
        charge = -hit->get_charge() * charge2nVs_; // we transform the charge to nV*s and take negative value to get positive charge
        assoc_gid = hit->get_geom_id();
      }
    }

    if(charge < charge_threshold_)
    {
      continue;
    }

    // vertex positions
    const geomtools::vector_3d & source_vertex_pos = source_vertex->get_spot().get_placement().get_translation();
    const geomtools::vector_3d & calo_vertex_pos = calo_vertex->get_spot().get_placement().get_translation();
    
    if(assoc_gid.get_type() == 1302) // set the "part" number for main wall modules
    {
      const geomtools::id_mgr & mgr = geo_manager_->get_id_mgr();
      mgr.set(assoc_gid, "part", 1);
    }

    int om_num = snemo::datamodel::om_num(assoc_gid);

    // calculate calo vertex position in the coordiante system of the OM
    const geomtools::mapping   & mapping = geo_manager_->get_mapping();
    const geomtools::geom_info & a_block_ginfo  = mapping.get_geom_info(assoc_gid);
    const geomtools::placement & a_block_world_placement = a_block_ginfo.get_world_placement();

	  geomtools::placement calo_vertex_pos_OM_placement;
    a_block_world_placement.relocate(calo_vertex_pos, calo_vertex_pos_OM_placement);
	  geomtools::vector_3d calo_vertex_pos_OM_vec = calo_vertex_pos_OM_placement.get_translation();

	  const double position_x =  (calo_vertex_pos_OM_vec.x());
    const double position_y =  (calo_vertex_pos_OM_vec.y());
    const double position_z = -(calo_vertex_pos_OM_vec.z());
    
    // save extracted values
    charge_ = charge;
    source_vertex_pos_.SetXYZ(source_vertex_pos.getX(), source_vertex_pos.getY(), source_vertex_pos.getZ());
    calo_vertex_pos_.SetXYZ(calo_vertex_pos.getX(), calo_vertex_pos.getY(), calo_vertex_pos.getZ());
    calo_vertex_pos_OM_.SetXYZ(position_y, -position_x, position_z);

    OM_data_[om_num]->Fill();
	}

  event_counter_++;

  return dpp::base_module::PROCESS_SUCCESS;
}

// check if the particle track has a vertex near a calibration source - if so, return the vertex
datatools::handle<snemo::datamodel::vertex> calibration_cuts_module::vertex_close_to_a_calib_source(const datatools::handle<snemo::datamodel::particle_track> & particle)
{
  for(const datatools::handle<snemo::datamodel::vertex> & vertex : particle->get_vertices())
  {
    if(vertex->is_on_reference_source_plane())
    {
      double ver_y = vertex->get_spot().get_position().getY();
      double ver_z = vertex->get_spot().get_position().getZ();
      // iterate through all calibration sources
      for(int i = 0;i < calib_source_rows_;i++)
      {
        for(int j = 0;j < calib_source_columns_;j++)
        {
          double ellipse_distance = ((ver_y - calib_source_Y_[i][j])*(ver_y - calib_source_Y_[i][j])) / (source_cut_ellipse_Y_*source_cut_ellipse_Y_) +
						((ver_z - calib_source_Z_[i][j])*(ver_z - calib_source_Z_[i][j])) / (source_cut_ellipse_Z_*source_cut_ellipse_Z_);

          if(ellipse_distance < 1.0) // vertex is close to a calibration source
          {
            return vertex;
          }
        }
      }
    }
  }
  return datatools::handle<snemo::datamodel::vertex>(nullptr);
}

// check if the particle track has exactly one associated OM hit - if so, return the vertex associated with the hit
datatools::handle<snemo::datamodel::vertex> calibration_cuts_module::track_has_one_assoc_calo(const datatools::handle<snemo::datamodel::particle_track> & particle)
{
  if(particle->get_associated_calorimeter_hits().size() == 1)
  {
    // iterate through vertices to find the one on OM
    for(const datatools::handle<snemo::datamodel::vertex> & vertex : particle->get_vertices())
    {
      // for now we ignore gveto
      if(vertex->is_on_main_calorimeter() || vertex->is_on_x_calorimeter())
      {
        return vertex;
      }
    }
  }
  return datatools::handle<snemo::datamodel::vertex>(nullptr);
}