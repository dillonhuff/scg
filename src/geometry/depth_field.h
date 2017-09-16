#pragma once

#include <cstdlib>

#include "geometry/point.h"
#include "utils/check.h"

namespace gca {

  class depth_field {
  protected:
    point origin;
    float* column_heights;

  public:
    double resolution;
    double x_len, y_len;
    int num_x_elems, num_y_elems;

    depth_field(const point p_origin,
		double x_w,
		double y_w,
		double xy_resolution) :
      origin(p_origin),
      resolution(xy_resolution), x_len(x_w), y_len(y_w),
      num_x_elems(x_w/static_cast<double>(xy_resolution)),
      num_y_elems(y_w/static_cast<double>(xy_resolution)) {

      DBG_ASSERT(origin.z == 0.0);

      column_heights =
	static_cast<float*>(malloc(sizeof(float)*num_x_elems*num_y_elems));

      for (int i = 0; i < num_x_elems; i++) {
	for (int j = 0; j < num_y_elems; j++) {
	  set_column_height(i, j, 0);
	}
      }
    }

    depth_field(const depth_field& other) :
      origin(other.origin),
      resolution(other.resolution), x_len(other.x_len), y_len(other.y_len),
      num_x_elems(other.num_x_elems),
      num_y_elems(other.num_y_elems) {

      DBG_ASSERT(origin.z == 0.0);

      column_heights =
	static_cast<float*>(malloc(sizeof(float)*num_x_elems*num_y_elems));

      for (int i = 0; i < num_x_elems; i++) {
	for (int j = 0; j < num_y_elems; j++) {
	  set_column_height(i, j, other.column_height(i, j));
	}
      }
    }
    
    inline point get_origin() const { return origin; }

    inline double x_center(const int i) const {
      return x_min() + resolution*i + (resolution/2.0);
    }

    inline double y_center(const int i) const {
      return y_min() + resolution*i + (resolution/2.0);
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

    double z_max() const {
      auto max_it =
	max_element(column_heights,
		    column_heights + num_x_elems*num_y_elems,
		    [](const float l, const float r) {
		      return l < r;
		    });

      DBG_ASSERT(max_it);

      return *max_it;
    }

    double z_min() const {
      auto min_it =
	min_element(column_heights,
		    column_heights + num_x_elems*num_y_elems,
		    [](const float l, const float r) {
		      return l < r;
		    });

      DBG_ASSERT(min_it);

      return *min_it;
    }
    
    inline double x_coord(const int i) const {
      return get_origin().x + resolution*i;
    }

    inline double y_coord(const int i) const {
      return get_origin().y + resolution*i;
    }
    
    inline int x_index(double x) const {
      return static_cast<int>((x - x_min()) / resolution);
      //return static_cast<int>(x / resolution);
    }

    inline int y_index(double y) const {
      return static_cast<int>((y - y_min()) / resolution);
    }
    
    void set_height(double x_s, double x_e,
		    double y_s, double y_e,
		    double h) {
      int first_x = x_index(x_s);
      int last_x = x_index(x_e);
      int first_y = y_index(y_s);
      int last_y = y_index(y_e);
      for (int i = first_x; i < last_x; i++) {
	for (int j = first_y; j < last_y; j++) {
	  set_column_height(i, j, h);
	}
      }      
    }
    
    ~depth_field() {
      delete[] column_heights;
    }

    inline float column_height(int i, int j) const {
      return *(column_heights + i*num_y_elems + j);
    }

    inline void set_column_height(int i, int j, float f) {
      *(column_heights + i*num_y_elems + j) = f;
    }

    inline bool legal_column(int i, int j) const {
      return (0 <= i && i < num_x_elems) && (0 <= j && j < num_y_elems);
    }

  };

}
