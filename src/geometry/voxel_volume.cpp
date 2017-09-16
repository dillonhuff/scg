#include "geometry/voxel_volume.h"

namespace gca {

  voxel_volume::voxel_volume(const point p_origin,
			     const double p_x_len,
			     const double p_y_len,
			     const double p_z_len,
			     const double p_resolution) :
    origin(p_origin),
    x_len(p_x_len),
    y_len(p_y_len),
    z_len(p_z_len),
    resolution(p_resolution),
    nx_elems(x_len / resolution),
    ny_elems(y_len / resolution),
    nz_elems(z_len / resolution) {
    
  }
  
}
