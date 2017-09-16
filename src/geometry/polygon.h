#ifndef GCA_POLYGON_H
#define GCA_POLYGON_H

#include "geometry/box.h"
#include "geometry/line.h"
#include "geometry/point.h"
#include "geometry/polyline.h"

namespace gca {

  class oriented_polygon {
  protected:
    vector<point> verts;
    point norm;

  public:

    oriented_polygon() :
      verts({}), norm(point(0, 0, 0)) {}

    // TODO: Add assertion about vertex values at start and end?
    oriented_polygon(point normalp, const vector<point>& verticesp) :
      verts(verticesp), norm(normalp) {
    }

    point normal() const { return norm; }

    inline const vector<point>& vertices() const { return verts; }

    inline point pt(unsigned i) const { return vertices()[i]; }
    inline double height() const
    { return vertices().front().z; }
  };

  bool contains(const oriented_polygon& g, point p);
  bool contains(const oriented_polygon& poly,
		const oriented_polygon& maybe_contained);

  bool is_horizontal(const oriented_polygon& p);

  bool contains(const oriented_polygon& poly,
		const oriented_polygon& maybe_contained);

  bool overlaps(line l, const oriented_polygon& p);

  box bounding_box(const oriented_polygon& p);

  oriented_polygon project(const oriented_polygon& p, double z);

  polyline to_polyline(const oriented_polygon& p);

  oriented_polygon extract_boundary(std::vector<oriented_polygon>& polygons);

  std::vector<polyline>
  clip_polyline_along(const polyline& p,
		      const std::vector<oriented_polygon>& holes);

  vector<polyline> sample_lines_2d(const oriented_polygon& b, double x_inc, double y_inc, double z_level);


  oriented_polygon base(const box b);

  double min_z(const oriented_polygon& p);

  double area(const oriented_polygon& p);
  double signed_area(const oriented_polygon& p);

}

#endif
