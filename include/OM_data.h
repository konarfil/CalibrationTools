#ifndef OM_DATA_H
#define OM_DATA_H

// Bayeux
#include <geomtools/geom_id.h>

// Other
#include <vector>

struct OM_data
{
    std::vector<double> charge; // charged measured by given OM
    std::vector<double> vertex_on_OM_X; // horizontal position of vertices relativelly to OM front face middle
    std::vector<double> vertex_on_OM_Y; // vertical position of vertices relativelly to OM front face middle
    std::vector<double> one_over_cos; // 1 over cosine of track angle
    geomtools::geom_id  gid; // GID of given OM
};
#endif