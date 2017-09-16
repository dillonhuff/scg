#ifndef GCA_BOX_H
#define GCA_BOX_H

#include <cassert>
#include <iostream>
#include <vector>

#include "geometry/point.h"
#include "geometry/polyline.h"

using namespace std;

namespace gca {

  struct box {
    double x_min, x_max, y_min, y_max, z_min, z_max;
    box(double x_minp, double x_maxp,
	double y_minp, double y_maxp,
	double z_minp, double z_maxp) :
      x_min(x_minp), x_max(x_maxp),
      y_min(y_minp), y_max(y_maxp),
      z_min(z_minp), z_max(z_maxp) {
      if (x_min > x_max)
	{ cout << "Bad X range: " << x_min << " > " << x_max << endl; assert(false); } 
      if (y_min > y_max)
	{ cout << "Bad Y range: " << y_min << " > " << y_max << endl; assert(false); } 
      if (z_min > z_max)
	{ cout << "Bad Z range: " << z_min << " > " << z_max << endl; assert(false); } 
    }

    bool contains(const point p) {
      return (x_min < p.x && p.x < x_max) &&
		      (y_min < p.y && p.y < y_max) &&
			       (z_min < p.z && p.z < z_max);
    }

    double x_len() const { return x_max - x_min; }
    double y_len() const { return y_max - y_min; }
    double z_len() const { return z_max - z_min; }
  };

  ostream& operator<<(ostream& out, const box& b);

  bool overlap(const box l, const box r);

  box bound_positions(const vector<point>& pts);
  box bound_boxes(const vector<box>& boxes);
  bool fits_inside(const box& inner, const box& outer);

  vector<point> sample_points_2d(const box b, double x_inc, double y_inc, double z_level);

  vector<polyline> sample_lines_2d(const box b, double x_inc, double y_inc, double z_level);
  vector<point> sample_points_3d(const box b,
				 double x_inc,
				 double y_inc,
				 double z_inc);

  template<typename F>
  vector<point> sample_filtered_points_2d(const box b,
					  double x_inc,
					  double y_inc,
					  double z_level,
					  F f) {
    auto points = sample_points_2d(b, x_inc, y_inc, z_level);
    delete_if(points, f);
    return points;
  }

  template<typename F>
  vector<point> sample_filtered_points_3d(const box b,
					  double x_inc,
					  double y_inc,
					  double z_inc,
					  F f) {
    auto points = sample_points_3d(b, x_inc, y_inc, z_inc);
    delete_if(points, f);
    return points;
  }
  
  template<typename InputIt>
  box bounding_box(InputIt s, InputIt e) {
    vector<box> boxes;
    while (s != e) {
      boxes.push_back(bounding_box(*s));
      ++s;
    }
    return bound_boxes(boxes);
  }

  // TODO: Add some unit tests
  template<typename T>
  bool intervals_overlap(pair<T, T> i, pair<T, T> j) {
    return (j.first <= i.second && i.second <= j.second) ||
      (j.first <= i.first && i.first <= j.second) ||
      (i.first <= j.first && i.first <= j.second &&
       j.first <= i.second && j.second <= i.second);
  }
  
}

#endif
