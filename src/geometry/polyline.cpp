#include <cmath>

#include "geometry/line.h"
#include "geometry/polyline.h"
#include "utils/algorithm.h"

namespace gca {

  bool is_closed(const polyline& p) {
    return within_eps(p.front(), p.back());
  }

  double signed_area(const polyline& p) {
    assert(is_closed(p));
    double signed_a = 0.0;
    for (auto l : p.lines()) {
      signed_a += (l.end.x - l.start.x)*(l.end.y + l.start.y);
    }
    return signed_a / 2.0;
  }

  double area(const polyline& p) {
    double sa = signed_area(p);
    return abs(sa);
  }

  offset_dir exterior_direction(const polyline& p) {
    assert(is_closed(p));
    return signed_area(p) >= 0.0 ? OFFSET_LEFT : OFFSET_RIGHT;
  }

  offset_dir interior_direction(const polyline& p) {
    return exterior_direction(p) == OFFSET_LEFT ? OFFSET_RIGHT : OFFSET_LEFT;
  }

  bool pointwise_within_eps(const polyline& p, const polyline& q, double tol) {
    if (!(p.num_points() == q.num_points())) { return false; }
    return mismatch(p.begin(), p.end(), q.begin(),
		    [tol](point l, point r) { return within_eps(l, r, tol); }).first == p.end();
  }

  // TODO: Add assertion to check for simplicity, e.g. the only
  // intersection of the lines is potentially the first and last point
  polyline offset(const polyline& p, offset_dir d, double n) {
    assert(p.num_points() > 1);
    vector<line> initial_lines;
    double degrees = d == OFFSET_LEFT ? 90.0 : -90.0;
    for (auto l : p.lines()) {
      point diff = l.end - l.start;
      diff.z = 0;
      auto offset_vec = n * (diff.normalize()).rotate_z(degrees);
      initial_lines.push_back(line(l.start + offset_vec, l.end + offset_vec));
    }
    vector<point> new_points(initial_lines.size());
    apply_between(initial_lines.begin(), initial_lines.end(),
		  new_points.begin() + 1,
		  trim_or_extend_unsafe);
    // TODO: Allow this tolerance to be set by the api
    if (!within_eps(p.front(), p.back(), 0.000001)) {
      new_points[0] = initial_lines[0].start;
      new_points.push_back(initial_lines.back().end);
    } else {
      assert(p.num_points() > 3);
      point n = trim_or_extend_unsafe(initial_lines[initial_lines.size() - 1],
				      initial_lines[0]);
      new_points[0] = n;
      new_points.push_back(n);
    }
    polyline off(new_points);
    return off;
  }

  vector<point> points(const polyline& p) {
    vector<point> pts;
    for (unsigned i = 0; i < p.num_points(); i++) {
      pts.push_back(p.pt(i));
    }
    return pts;
  }

  vector<point> points(const vector<polyline>& p) {
    vector<point> pts;
    for (auto pl : p) {
      auto pl_pts = points(pl);
      pts.insert(end(pts), begin(pl_pts), end(pl_pts));
    }
    return pts;
  }

  polyline shift(const polyline& p, const point s) {
    std::vector<point> pts;
    for (auto pt : p) {
      pts.push_back(pt + s);
    }
    return polyline(pts);
  }
  
  std::vector<polyline> shift_lines(const std::vector<polyline>& lines,
				    const point s) {
    vector<polyline> slines;
    for (auto pl : lines) {
      slines.push_back(shift(pl, s));
    }
    return slines;
  }

  double min_in_dir(const std::vector<polyline>& lines, const point dir) {
    return min_distance_along(points(lines), dir);
  }

  double max_in_dir(const std::vector<polyline>& lines, const point dir) {
    assert(lines.size() > 0);
    return max_distance_along(points(lines), dir);
  }

  polyline project_line(const polyline& l, const double d) {
    vector<point> pts;
    for (auto p : l) {
      pts.push_back(point(p.x, p.y, d));
    }
    return pts;
  }

  std::vector<polyline> project_lines(const std::vector<polyline>& lines,
				      const double d) {
    vector<polyline> ls;
    for (auto line : lines) {
      ls.push_back(project_line(line, d));
    }
    return ls;
  }

  double length(const polyline& pl) {
    DBG_ASSERT(pl.num_points() > 1);

    double length = 0.0;

    for (unsigned i = 0; i < pl.num_points() - 1; i++) {
      length += (pl.pt(i + 1) - pl.pt(i)).len();
    }

    return length;
  }
}
