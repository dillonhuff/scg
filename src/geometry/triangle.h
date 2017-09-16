#ifndef GCA_TRIANGLE_H
#define GCA_TRIANGLE_H

#include <cmath>
#include <vector>

#include "geometry/box.h"
#include "geometry/point.h"
#include "geometry/polygon.h"

using namespace std;

namespace gca {

  struct triangle {
    point normal;
    point v1;
    point v2;
    point v3;

    triangle(point normalp, point v1p, point v2p, point v3p) :
      normal(normalp), v1(v1p), v2(v2p), v3(v3p) {}

    std::vector<line> edges() const {
      return {line(v1, v2), line(v2, v3), line(v3, v1)};
    }

    inline point centroid() const
    { return (1.0 / 3.0) * (v1 + v2 + v3); }

    double area() const {
      double a = (v1  - v2).len();
      double b = (v2  - v3).len();
      double c = (v1  - v3).len();
      double s = (a + b + c) / 2.0;
      return sqrt(s*(s - a)*(s - b)*(s - c));
    }
  };

  bool in_projection(const triangle t, const point p);
  bool below(const triangle t, const point p);
  double z_at(const triangle t, double x, double y);
  bool intersects(const triangle t, const line l);

  double min_z(const std::vector<triangle>& triangles);
  bool is_upward_facing(const triangle& t, double tolerance);
  bool same_orientation(const triangle& x, const triangle& y, double tolerance);

  // TODO: Change this to operate on a triangular mesh
  std::vector<oriented_polygon> mesh_bounds(const std::vector<triangle>& tris);

  vector<oriented_polygon>
  unordered_segments_to_polygons(point normal,
				 vector<line>& lines);
  
  ostream& operator<<(ostream& out, const triangle& t);

  bool intersects_triangles(line l, const std::vector<triangle>& triangles);

  double distance_along(point normal, const triangle t);

  std::vector<triangle>
  triangulate_box_pts(const std::vector<point>& p);

  std::vector<triangle> box_triangles(box b);


  std::vector<triangle> triangulate_polygon(const oriented_polygon& p);

  double max_z(const triangle t);

  point normal(const triangle t);

}

#endif
