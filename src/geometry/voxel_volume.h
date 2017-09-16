#pragma once

#include <boost/container/flat_set.hpp>

#include "geometry/point.h"

namespace gca {

  struct voxel {
    int xi, yi, zi;
  };

  static inline bool operator<(const voxel l, const voxel r) {
    if (l.xi < r.xi) { return true; }
    if (l.xi > r.xi) { return false; }

    if (l.yi < r.yi) { return true; }
    if (l.yi > r.yi) { return false; }

    if (l.zi < r.zi) { return true; }

    return false;
      
  }

  class voxel_volume {

  protected:

    point origin;
    double x_len, y_len, z_len;
    double resolution;
    int nx_elems, ny_elems, nz_elems;

    boost::container::flat_set<voxel> voxels;
    
  public:
    voxel_volume(const point origin,
		 const double x_len,
		 const double y_len,
		 const double z_len,
		 const double resolution);

    inline double get_resolution() const { return resolution; }

    inline point get_origin() const { return origin; }

    inline bool is_empty(int x_i, int y_i, int z_i) const {
      return !is_occupied(x_i, y_i, z_i);
    }

    inline bool is_occupied(int x_i, int y_i, int z_i) const {
      return voxels.find( voxel{x_i, y_i, z_i} ) != voxels.end();
    }

    inline void set_occupied(int x_i, int y_i, int z_i) {
      voxels.insert( voxel{x_i, y_i, z_i} );
    }

    inline double x_center(const int i) const {
      return x_min() + resolution*i + (resolution/2.0);
    }

    inline double y_center(const int i) const {
      return y_min() + resolution*i + (resolution/2.0);
    }

    inline double z_center(const int i) const {
      return z_min() + resolution*i + (resolution/2.0);
    }

    inline double x_length() const {
      return x_max() - x_min();
    }

    inline double y_length() const {
      return y_max() - y_min();
    }

    inline double z_length() const {
      return z_max() - z_min();
    }
    
    inline double x_min() const {
      return origin.x;
    }

    inline double x_max() const {
      return origin.x + x_len;
    }

    inline double y_min() const {
      return origin.y;
    }

    inline double y_max() const {
      return origin.y + y_len;
    }

    inline double z_min() const {
      return origin.z;
    }
    
    inline double z_max() const {
      return origin.z + z_len;
    }

    inline int num_x_elems() const {
      return nx_elems;
    }

    inline int num_y_elems() const {
      return ny_elems;
    }

    inline int num_z_elems() const {
      return nz_elems;
    }
    
  };


}
