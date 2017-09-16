#include <cmath>

#include "geometry/polygon.h"
#include "geometry/polyline.h"
#include "geometry/vtk_debug.h"
#include "utils/algorithm.h"

namespace gca {

  /*
    Return the angle between two vectors on a plane
    The angle is from vector 1 to vector 2, positive anticlockwise
    The result is between -pi -> pi
  */
  double angle_2d(double x1, double y1, double x2, double y2)
  {
    double dtheta, theta1, theta2;

    theta1 = atan2(y1,x1);
    theta2 = atan2(y2,x2);
    dtheta = theta2 - theta1;
    while (dtheta > M_PI)
      dtheta -= 2*M_PI; //TWOPI;
    while (dtheta < -M_PI)
      dtheta += 2*M_PI; //TWOPI;

    return(dtheta);
  }

  bool contains(const oriented_polygon& poly,
		const oriented_polygon& maybe_contained) {
    for (auto pt : maybe_contained.vertices()) {
      if (!contains(poly, pt)) {
	return false;
      }
    }
    return true;
  }

  bool contains(const oriented_polygon& poly, point p)
  {
    double angle = 0;
    point p1, p2;
    int n = poly.vertices().size();

    for (int i = 0; i < n; i++) {
      p1.x = poly.pt(i).x - p.x;
      p1.y = poly.pt(i).y - p.y;
      p2.x = poly.pt((i+1)%n).x - p.x;
      p2.y = poly.pt((i+1)%n).y - p.y;
      angle += angle_2d(p1.x, p1.y, p2.x, p2.y);
    }

    if (abs(angle) < M_PI)
      { return false; }
    else
      { return true; }
  }

  bool is_horizontal(const oriented_polygon& p) {
    return within_eps(p.normal().z, 1.0, 0.001);
  }

  box bounding_box(const oriented_polygon& p) {
    return bound_positions(p.vertices());
  }

  bool overlaps(line l, const oriented_polygon& p) {
    polyline pl(p.vertices());
    for (auto pll : pl.lines()) {
      if (segment_intersection_2d(l, pll).just)
	{ return true; }
    }
    return false;
  }

  template<typename InputIt>
  bool overlaps_any(line l, InputIt s, InputIt e) {
    while (s != e) {
      if (overlaps(l, *s)) {
	return true;
      }
      ++s;
    }
    return false;
  }

  oriented_polygon project(const oriented_polygon& p, double z) {
    vector<point> pts;
    for (auto pt : p.vertices()) {
      pts.push_back(point(pt.x, pt.y, z));
    }
    return oriented_polygon(p.normal(), pts);
  }

  polyline to_polyline(const oriented_polygon& p) {
    auto v = p.vertices();
    v.push_back(p.vertices().front());
    return polyline(v);
  }

  oriented_polygon extract_boundary(std::vector<oriented_polygon>& polygons) {
    DBG_ASSERT(polygons.size() > 0);

    for (unsigned i = 0; i < polygons.size(); i++) {
      auto possible_bound = polygons[i];
      bool contains_all = true;
      for (unsigned j = 0; j < polygons.size(); j++) {
	if (i != j) {
	  auto possible_hole = polygons[j];
	  if (!contains(possible_bound, possible_hole)) {
	    contains_all = false;
	    break;
	  }
	}
      }
      if (contains_all) {
	polygons.erase(polygons.begin() + i);
	return possible_bound;
      }
    }

    vtk_debug_polygons(polygons);

    DBG_ASSERT(false);
  }

  vector<polyline> clip_polyline_along(const polyline& p,
				       const vector<oriented_polygon>& holes) {
    auto inside_any =
      [holes](const point l) { 
      return any_of(begin(holes), end(holes),
		    [l](const oriented_polygon& pg) {
		      return contains(pg, l);
		    });
    };

    auto intersects_none =
      [holes](const point a, const point b) {
      return !any_of(begin(holes), end(holes),
		    [a, b](const oriented_polygon& pg) {
		      return overlaps(line(a, b), pg);
		     });
    };

    vector<point> pts(begin(p), end(p));
    // TODO: Deal with convex case
    delete_if(pts, inside_any);
    vector<vector<point>> lgs = split_by(pts, intersects_none);
    vector<polyline> lines;
    for (auto l : lgs) {
      lines.push_back(polyline(l));
    }
    return lines;
  }

  vector<polyline> sample_lines_2d(const oriented_polygon& b, double x_inc, double y_inc, double z_level) {
    vector<polyline> box_lines =
      sample_lines_2d(bounding_box(b), x_inc, y_inc, z_level);
    vector<polyline> lines;
    for (auto bl : box_lines) {
      vector<point> pts(begin(bl), end(bl));
      delete_if(pts,
		[b](const point p)
		{ return !contains(b, p); });
      lines.push_back(pts);
    }
    return lines;
  }

  oriented_polygon base(const box b) {
    point p1(b.x_min, b.y_min, b.z_min);
    point p2(b.x_min, b.y_max, b.z_min);
    point p3(b.x_max, b.y_max, b.z_min);
    point p4(b.x_max, b.y_min, b.z_min);
    vector<point> verts{p1, p2, p3, p4};
    oriented_polygon interior(point(0, 0, 1), verts);
    return interior;
  }

  double min_z(const oriented_polygon& p) {
    auto min_z =
      max_element(begin(p.vertices()), end(p.vertices()),
		  [](const point l, const point r)
		  { return l.z < r.z; });
    return (*min_z).z;
  }

  double area(const oriented_polygon& p) {
    return area(to_polyline(p));
  }
  
  double signed_area(const oriented_polygon& p) {
    return signed_area(to_polyline(p));
  }

}
