#pragma once

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include "geometry/box.h"
#include "geometry/plane.h"
#include "geometry/ring.h"

namespace gca {

  namespace bg = boost::geometry;

  typedef bg::model::d2::point_xy<double> boost_point_2;
  typedef bg::model::polygon<boost_point_2> boost_poly_2;
  typedef bg::model::multi_polygon<boost_poly_2> boost_multipoly_2;
  typedef bg::model::multi_point<boost_point_2> boost_multipoint_2;

  typedef bg::model::linestring<boost_point_2> boost_linestring_2;
  typedef bg::model::multi_linestring<boost_linestring_2> boost_multilinestring_2;

  class polygon_3 {
  protected:
    std::vector<point> outer_ring;
    std::vector<std::vector<point>> inner_rings;
  public:
    // polygon_3(const std::vector<point> vertices) :
    //   outer_ring(vertices) {

    //   outer_ring = clean_vertices(outer_ring);
    //   delete_antennas(outer_ring);

    //   DBG_ASSERT(outer_ring.size() >= 3);
    // }

    polygon_3(const std::vector<point>& vertices,
	      const std::vector<std::vector<point>>& hole_verts,
	      const bool dummy_param);
    
    polygon_3(const std::vector<point>& vertices, const bool dummy_param) :
      outer_ring(vertices) {

      DBG_ASSERT(outer_ring.size() >= 3);
    }
    
    // polygon_3(const std::vector<point> vertices,
    // 		      const std::vector<std::vector<point>> hole_verts);

    void add_hole(const std::vector<point>& hole) {
      inner_rings.push_back(hole);
    }

    point normal() const {
      return ring_normal(outer_ring);
    }

    point vertex(const unsigned i) const { return outer_ring[i]; }

    const std::vector<point>& vertices() const { return outer_ring; }

    const std::vector<std::vector<point>>& holes() const { return inner_rings; }

    std::vector<point> hole(unsigned i) const { return inner_rings[i]; }
    
    unsigned num_vertices() const { return outer_ring.size(); }

    void correct_winding_order(const point n) {
      gca::correct_winding_order(outer_ring, n);
      for (std::vector<point>& ir : inner_rings) {
	gca::correct_winding_order(ir, n);
      }
    }
  };

  typedef polygon_3 labeled_polygon_3;

  polygon_3 shift(const point p, const labeled_polygon_3& poly);
  std::vector<polygon_3> shift(const point p,
			       const std::vector<polygon_3>& poly);

  void check_simplicity(const labeled_polygon_3& p);

  typedef std::vector<std::vector<labeled_polygon_3>> surface_levels;

  polygon_3 smooth_buffer(const polygon_3& p,
			  const double tol);

  labeled_polygon_3 dilate(const labeled_polygon_3& p, const double tol);

  boost_poly_2
  to_boost_poly_2(const labeled_polygon_3& p);

  labeled_polygon_3 project_onto(const plane p,
				 const labeled_polygon_3& poly);

  std::vector<point> project_points(const plane pl,
				    const std::vector<point>& pts);

  boost_poly_2 rotate_to_2D(const labeled_polygon_3& p);

  double area(const polygon_3& p);

  double area(const std::vector<polygon_3>& p);

  std::vector<polygon_3>
  planar_polygon_union(const std::vector<polygon_3>& polys);

  std::vector<polygon_3>
  planar_polygon_union(const point n,
		       const std::vector<polygon_3>& polys);
  
  boost_multipoly_2
  planar_union_boost(const std::vector<polygon_3>& polys);
  
  polygon_3
  convex_hull_2D(const std::vector<point>& pts,
		 const point n,
		 const double z_level);

  polygon_3
  to_polygon_3(const double z, const boost_poly_2& p);

  box bounding_box(const polygon_3& p);

  polygon_3 project(const polygon_3& p, double z);

  boost_multipoly_2
  to_boost_multipoly_2(const std::vector<polygon_3>& lines);

  polygon_3
  build_clean_polygon_3(const std::vector<point>& vertices,
			const std::vector<std::vector<point>>& hole_verts);

  polygon_3
  build_clean_polygon_3(const std::vector<point>& vertices);

  polyline to_polyline(const polygon_3& p);

  std::vector<polyline> to_polylines(const polygon_3& p);

  std::vector<polygon_3>
  arrange_rings(const std::vector<std::vector<point>>& rings);

  labeled_polygon_3 shrink(const labeled_polygon_3& p, const double tol);

  boost::optional<labeled_polygon_3>
  shrink_optional(const labeled_polygon_3& p,
		  const double tol);

  std::vector<polygon_3> shift(const point p,
			       const std::vector<polygon_3>& polys);


  std::vector<labeled_polygon_3>
  dilate_polygons(const std::vector<labeled_polygon_3>& polys, const double tol);

  int curve_count(const polygon_3& f);

  // vector<polygon_3> from_boost_multipoly_2(const boost_multipoly_2& p,
  // 					   const rotation& r,
  // 					   const double z_level);
  
  std::vector<polygon_3>
  polygon_difference(const std::vector<polygon_3>& as,
		     const std::vector<polygon_3>& bs);

  bool
  contains(const std::vector<polygon_3>& as,
	   const std::vector<polygon_3>& bs);

  std::vector<polygon_3>
  polygon_union(const std::vector<polygon_3>& as,
		const std::vector<polygon_3>& bs);
  
  std::vector<polygon_3>
  polygon_intersection(const std::vector<polygon_3>& as,
		       const std::vector<polygon_3>& bs);
  
}
