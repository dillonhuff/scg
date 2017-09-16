#ifndef GCA_POLYLINE_H
#define GCA_POLYLINE_H

#include <cassert>
#include <vector>

#include "geometry/line.h"
#include "geometry/point.h"

using namespace std;

namespace gca {

  class polyline {
  protected:
    vector<point> points;

  public:

    polyline(const vector<point>& pointsp) : points(pointsp) {}

    unsigned num_points() const { return points.size(); }
    vector<point>::iterator begin() { return points.begin(); }
    vector<point>::iterator end() { return points.end(); }
    vector<point>::const_iterator begin() const { return points.begin(); }
    vector<point>::const_iterator end() const { return points.end(); }

    point front() const { return points.front(); }
    point back() const { return points.back(); }

    point pt(unsigned n) const { return points[n]; }

    vector<line> lines() const {
      assert(points.size() > 1);
      vector<line> ls;
      for (vector<point>::const_iterator it = points.begin();
	   it != points.end() - 1; ++it)
	{ ls.push_back(line(*it, *(it + 1))); }
      return ls;
    }

    template<typename F>
    polyline apply(F f) const {
      std::vector<point> pts = points;
      std::transform(std::begin(points), std::end(points), std::begin(pts), f);
      return polyline(pts);
    }
      
  };

  template<typename F>
  polyline apply_to_points(const polyline& l, F f) {
    vector<point> pts;
    for (auto p : l) {
      pts.push_back(f(p));
    }
    return polyline(pts);
  }

  enum offset_dir {
    OFFSET_LEFT,
    OFFSET_RIGHT
  };

  bool is_closed(const polyline& p);
  offset_dir exterior_direction(const polyline& p);
  offset_dir interior_direction(const polyline& p);
  bool pointwise_within_eps(const polyline& p, const polyline& q, double tol);

  double area(const polyline& p);
  double signed_area(const polyline& p);

  polyline offset(const polyline& p, offset_dir d, double n);

  std::vector<point> points(const polyline& p);
  std::vector<point> points(const vector<polyline>& p);

  std::vector<polyline> shift_lines(const std::vector<polyline>& lines,
				    const point s);

  polyline shift(const polyline& p, const point s);

  double max_in_dir(const std::vector<polyline>& lines, const point dir);
  double min_in_dir(const std::vector<polyline>& lines, const point dir);

  std::vector<polyline> project_lines(const std::vector<polyline>& lines,
				      const double d);

  double length(const polyline& pl);
}

#endif
