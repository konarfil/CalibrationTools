#ifndef CALIBRATION_CUTS_DRIVER_H
#define CALIBRATION_CUTS_DRIVER_H

#include "OM_data.h"

// Falaise
#include <falaise/snemo/services/service_handle.h>
#include <falaise/snemo/services/geometry.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/datamodels/particle_track.h>
#include <falaise/snemo/datamodels/precalibrated_data.h>

// Bayeux
#include <geomtools/manager.h>
#include <geomtools/geometry_service.h>

// Other
#include <map>
#include <memory>

class calibration_cuts_driver
{
  public:
    // apply data cuts - if the event passes, save its information into OM data
    void process_event(const snemo::datamodel::particle_track_data & PTD, const snemo::datamodel::precalibrated_data & pCD);

    // Constructor 
    calibration_cuts_driver(std::shared_ptr<std::map<int, OM_data>> OM_data, const snemo::service_handle<snemo::geometry_svc> & geo_manager, 
                            double source_cut_ellipse_y, double source_cut_ellipse_z);

    // Destructor
    ~calibration_cuts_driver();

  private:
    // check if the particle track has a vertex near a calibration source - if so, return the vertex
    datatools::handle<snemo::datamodel::vertex> vertex_close_to_a_calib_source(const datatools::handle<snemo::datamodel::particle_track> & particle);

    // check if the particle track has exactly one associated OM hit - if so, return the vertex associated with the hit
    datatools::handle<snemo::datamodel::vertex> track_has_one_assoc_calo(const datatools::handle<snemo::datamodel::particle_track> & particle);

    const snemo::service_handle<snemo::geometry_svc> & geo_manager_;

    std::shared_ptr<std::map<int, OM_data>> OM_data_; // data corresponding to indivitual OMs

    static constexpr double charge2nVs = 1e6; // multiplicative constant to transform charge to units nV*s

    // number of calibration sources
    static const int calib_source_rows_ = 7;
    static const int calib_source_columns_ = 6;

    // dimensions of the ellipse for calibration source vertex cut
    double source_cut_ellipse_Y_;
    double source_cut_ellipse_Z_;

    // Y and Z coorinates of centres of calibration sources
    double calib_source_Y_[calib_source_rows_][calib_source_columns_];
    double calib_source_Z_[calib_source_rows_][calib_source_columns_];
};
#endif