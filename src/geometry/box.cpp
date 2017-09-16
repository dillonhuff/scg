#include <iostream>

#include "geometry/box.h"

using namespace std;

namespace gca {

  ostream& operator<<(ostream& out, const box& b) {
    out << "X min = " << b.x_min << endl;
    out << "X max = " << b.x_max << endl;
    out << "Y min = " << b.y_min << endl;
    out << "Y max = " << b.y_max << endl;
    out << "Z min = " << b.z_min << endl;
    out << "Z max = " << b.z_max << endl;
    return out;
  }

  bool overlap(const box l, const box r) {
    bool x_overlap = intervals_overlap(pair<double, double>(l.x_min, l.x_max),
				       pair<double, double>(r.x_min, r.x_max));
    bool y_overlap = intervals_overlap(pair<double, double>(l.y_min, l.y_max),
				       pair<double, double>(r.y_min, r.y_max));
    bool z_overlap = intervals_overlap(pair<double, double>(l.z_min, l.z_max),
				       pair<double, double>(r.z_min, r.z_max));
    return x_overlap && y_overlap && z_overlap;
  }

  box bound_positions(const vector<point>& pts) {
    assert(pts.size() > 0);
    auto xminmax = minmax_element(pts.begin(), pts.end(),
				  [](const point l, const point r)
				  { return l.x < r.x; });
    assert(xminmax.first != pts.end());
    auto yminmax = minmax_element(pts.begin(), pts.end(),
				  [](const point l, const point r)
				  { return l.y < r.y; });
    auto zminmax = minmax_element(pts.begin(), pts.end(),
				  [](const point l, const point r)
				  { return l.z < r.z; });
    assert(yminmax.first != pts.end());
    box b((*(xminmax.first)).x, (*(xminmax.second)).x,
	  (*(yminmax.first)).y, (*(yminmax.second)).y,
	  (*(zminmax.first)).z, (*(zminmax.second)).z);
    return b;
  }

  box bound_boxes(const vector<box>& boxes) {
    //assert(boxes.size() > 0);
    if (boxes.size() == 0) {
      return box(0, 0, 0, 0, 0, 0);
    }

    auto x_min = *min_element(boxes.begin(), boxes.end(),
			      [](const box& l, const box& r)
			      { return l.x_min < r.x_min; });
    auto y_min = *min_element(boxes.begin(), boxes.end(),
			      [](const box& l, const box& r)
			      { return l.y_min < r.y_min; });
    auto z_min = *min_element(boxes.begin(), boxes.end(),
			      [](const box& l, const box& r)
			      { return l.z_min < r.z_min; });
    auto x_max = *max_element(boxes.begin(), boxes.end(),
			      [](const box& l, const box& r)
			      { return l.x_max < r.x_max; });
    auto y_max = *max_element(boxes.begin(), boxes.end(),
			      [](const box& l, const box& r)
			      { return l.y_max < r.y_max; });
    auto z_max = *max_element(boxes.begin(), boxes.end(),
			      [](const box& l, const box& r)
			      { return l.z_max < r.z_max; });
    return box(x_min.x_min, x_max.x_max,
	       y_min.y_min, y_max.y_max,
	       z_min.z_min, z_max.z_max);
  }

  bool fits_inside(const box& inner, const box& outer) {
    double inner_x_range = inner.x_max - inner.x_min;
    double inner_y_range = inner.y_max - inner.y_min;
    double inner_z_range = inner.z_max - inner.z_min;

    double outer_x_range = outer.x_max - outer.x_min;
    double outer_y_range = outer.y_max - outer.y_min;
    double outer_z_range = outer.z_max - outer.z_min;

    return (inner_x_range < outer_x_range) &&
      (inner_y_range < outer_y_range) &&
      (inner_z_range < outer_z_range);
  }

  vector<point> sample_points_2d(const box b, double x_inc, double y_inc, double z_level) {
    vector<point> pts;
    double x = b.x_min;
    while (x < b.x_max) {
      double y = b.y_min;
      while (y < b.y_max) {
	pts.push_back(point(x, y, z_level));
	y += y_inc;
      }
      x += x_inc;
    }
    return pts;
  }

  vector<polyline> sample_lines_2d(const box b, double x_inc, double y_inc, double z_level) {
    vector<polyline> pts;
    double x = b.x_min;
    while (x < b.x_max) {
      double y = b.y_min;
      vector<point> pline;
      while (y < b.y_max) {
	pline.push_back(point(x, y, z_level));
	y += y_inc;
      }
      pts.push_back(pline);
      x += x_inc;
    }
    return pts;
  }
  
  vector<point> sample_points_3d(const box b,
				 double x_inc,
				 double y_inc,
				 double z_inc) {
    vector<point> pts;
    double x = b.x_min;
    while (x < b.x_max) {
      double y = b.y_min;
      while (y < b.y_max) {
	double z = b.z_min;
	while (z < b.z_max) {
	  pts.push_back(point(x, y, z));
	  z += z_inc;
	}
	y += y_inc;
      }
      x += x_inc;
    }
    return pts;
  }

}
