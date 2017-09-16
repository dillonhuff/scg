#include "geometry/spline_sampling.h"

namespace gca {
  
  void append_spline(const b_spline* s,
		     vector<polyline>& polys) {
    // TODO: Come up with a better sampling strategy
    unsigned points_per_spline = 100;
    double next = 0.0;
    double inc = 1.0 / points_per_spline;
    vector<point> pts;
    for (unsigned i = 1; i < points_per_spline; i++) {
      next = next + inc;
      pts.push_back(s->eval(next));
    }
    polys.push_back(polyline(pts));
  }
  
  void append_splines(const vector<b_spline*>& splines,
		      vector<polyline>& polys)
  { for (auto s : splines) { append_spline(s, polys); } }
}
