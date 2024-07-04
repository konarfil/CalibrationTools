#include "calibration_cuts_driver.h"

// Falaise
#include <falaise/snemo/datamodels/geomid_utils.h>

// Bayeux
#include <geomtools/geometry_service.h>

// Other
#include <cmath>

// Constructor
calibration_cuts_driver::calibration_cuts_driver(std::shared_ptr<std::map<int, OM_data>> OM_data, const snemo::service_handle<snemo::geometry_svc> & geo_manager, 
                                                 double source_cut_ellipse_Y, double source_cut_ellipse_Z)
: OM_data_(OM_data), geo_manager_(geo_manager)
{
  source_cut_ellipse_Y_ = source_cut_ellipse_Y;
  source_cut_ellipse_Z_ = source_cut_ellipse_Z;

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
      const geomtools::vector_3d & calib_spot_pos  = calib_spot_placement.get_translation();

      calib_source_Y_[i][j] = calib_spot_pos.getY();
      calib_source_Z_[i][j] = calib_spot_pos.getZ();
    }
  }
}

// apply data cuts - if the event passes, save its information into OM data
void calibration_cuts_driver::process_event(const snemo::datamodel::particle_track_data & PTD, const snemo::datamodel::precalibrated_data & pCD)
{
  // iterate through all particles and find those coming from a calibration source and hitting an OM
  for (const datatools::handle<snemo::datamodel::particle_track> & particle : PTD.particles())
	{
    datatools::handle<snemo::datamodel::vertex> first_vertex = vertex_close_to_a_calib_source(particle);
    if(!first_vertex.has_data()) // there is no vertex close to a calibration source
    {
      continue;
    }

    datatools::handle<snemo::datamodel::vertex> second_vertex = track_has_one_assoc_calo(particle);
    if(!second_vertex.has_data()) // track does not have exactly one associated OM hit
    {
      continue;
    }

    // find the associated OM in pCD and save measured charge and its gid
    double charge = -1.0;
    geomtools::geom_id assoc_gid;
    for(const datatools::handle<snemo::datamodel::precalibrated_calorimeter_hit> & hit : pCD.calorimeter_hits())
    {
      if(hit->get_geom_id() == second_vertex->get_geom_id())
      {
        charge = -hit->get_charge() * charge2nVs; // we transform the charge to nV*s and take negative value to get positive charge
        assoc_gid = hit->get_geom_id();
      }
    }

    if(charge == -1.0) // probably should never happen
    {
      continue;
    }

    // vertex positions
    const geomtools::vector_3d & first_pos = first_vertex->get_spot().get_placement().get_translation();
    const geomtools::vector_3d & second_pos = second_vertex->get_spot().get_placement().get_translation();
    
    if(assoc_gid.get_type() == 1302) // set the "part" number for main wall modules
    {
      const geomtools::id_mgr & mgr = geo_manager_->get_id_mgr();
      mgr.set(assoc_gid, "part", 1);
    }

    // find associated OM position
    const geomtools::mapping & mapping = geo_manager_->get_mapping();
    const geomtools::geom_info & assoc_geom_info = mapping.get_geom_info(assoc_gid);
    const geomtools::placement & assoc_placement = assoc_geom_info.get_world_placement();
    const geomtools::vector_3d & assoc_pos  = assoc_placement.get_translation();

    //calculate horizontal and vertical distance between vertex and the center of associated OM
    double vertex_on_OM_X = -200.0;
    double vertex_on_OM_Y = -200.0;
    if(assoc_gid.get_type() == 1302)
    {
      vertex_on_OM_X = second_pos.getY() - assoc_pos.getY();
      vertex_on_OM_Y = second_pos.getZ() - assoc_pos.getZ();
    }
    else if(assoc_gid.get_type() == 1232)
    {
      vertex_on_OM_X = second_pos.getX() - assoc_pos.getX();
      vertex_on_OM_Y = second_pos.getZ() - assoc_pos.getZ();
    }
    else
    {
      continue;
    }

    // 1 over cosine of track angle = sqrt(1 + ((y_2 - y_1) / (x_2 - x_1))^2 + ((z_2 - z_1) / (x_2 - x_1))^2)
    // x1, y1, z1 - coordinates of the first vertex, x2, y2, z2 - coordinates of the second vertex
    double one_over_cos = std::sqrt(1.0 + (second_pos.getY() - first_pos.getY()) / ((second_pos.getX() - first_pos.getX())) * 
                                          (second_pos.getY() - first_pos.getY()) / ((second_pos.getX() - first_pos.getX())) + 
                                          (second_pos.getZ() - first_pos.getZ()) / ((second_pos.getX() - first_pos.getX())) * 
                                          (second_pos.getZ() - first_pos.getZ()) / ((second_pos.getX() - first_pos.getX())));

    // save information about the track into OM_data_
    int om_num = snemo::datamodel::om_num(assoc_gid);
    if((*OM_data_).find(om_num) == (*OM_data_).end())
    {
      (*OM_data_)[om_num] = OM_data{{}, {}, {}, {}, assoc_gid};
    }
    (*OM_data_)[om_num].charge.push_back(charge);
    (*OM_data_)[om_num].vertex_on_OM_X.push_back(vertex_on_OM_X);
    (*OM_data_)[om_num].vertex_on_OM_Y.push_back(vertex_on_OM_Y);
    (*OM_data_)[om_num].one_over_cos.push_back(one_over_cos);
	}
}

// check if the particle track has a vertex near a calibration source - if so, return the vertex
datatools::handle<snemo::datamodel::vertex> calibration_cuts_driver::vertex_close_to_a_calib_source(const datatools::handle<snemo::datamodel::particle_track> & particle)
{
  for(const datatools::handle<snemo::datamodel::vertex> & vertex : particle->get_vertices())
  {
    if(vertex->is_on_reference_source_plane())
    {
      double ver_y = vertex->get_spot().get_position().getY();
      double ver_z = vertex->get_spot().get_position().getZ();
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
datatools::handle<snemo::datamodel::vertex> calibration_cuts_driver::track_has_one_assoc_calo(const datatools::handle<snemo::datamodel::particle_track> & particle)
{
  if(particle->get_associated_calorimeter_hits().size() == 1)
  {
    // iterate through vertices to find the one on OM
    for(const datatools::handle<snemo::datamodel::vertex> & vertex : particle->get_vertices())
    {
      if(vertex->is_on_main_calorimeter() || vertex->is_on_x_calorimeter())
      {
        return vertex;
      }
    }
  }
  return datatools::handle<snemo::datamodel::vertex>(nullptr);
}