#include "geometry/polygon_3.h"

#include "geometry/rotation.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>


namespace gca {

  namespace bg = boost::geometry;

  typedef bg::model::d2::point_xy<double> boost_point_2;
  typedef bg::model::polygon<boost_point_2> boost_poly_2;
  //typedef bg::model::multi_polygon<boost_poly_2> boost_multipoly_2;
  //typedef bg::model::multi_point<boost_point_2> boost_multipoint_2;

  typedef bg::model::linestring<boost_point_2> boost_linestring_2;
  //typedef bg::model::multi_linestring<boost_linestring_2> boost_multilinestring_2;

  polygon_3
  build_clean_polygon_3(const std::vector<point>& vertices,
			const std::vector<std::vector<point>>& hole_verts) {
    vector<point> outer_ring = clean_vertices(vertices);
    delete_antennas(outer_ring);

    // There is an occasional test failure here in simple box
    if (!(outer_ring.size() >= 3)) {
      cout << "ERROR: Outer ring size = " << outer_ring.size() << endl;
      //vtk_debug_ring(outer_ring);

      DBG_ASSERT(outer_ring.size() >= 3);
    }

    vector<vector<point>> inner_rings;
    for (auto h : hole_verts) {

      auto new_h = clean_vertices(h);
      delete_antennas(new_h);

      if (!(new_h.size() >= 3)) {
	cout << "ERROR: Inner ring size = " << h.size() << endl;
      } else {
	inner_rings.push_back(new_h);
      }
    }

    return polygon_3(outer_ring, inner_rings, true);
  }

  polygon_3
  build_clean_polygon_3(const std::vector<point>& vertices) {
    return build_clean_polygon_3(vertices, {});
  }


//   // NOTE: These should be references
  polygon_3::polygon_3(const std::vector<point>& vertices,
		       const std::vector<std::vector<point>>& hole_verts,
		       const bool dummy_param) :
    outer_ring(vertices),
    inner_rings{} {

    for (auto h : hole_verts) {
      inner_rings.push_back(h);
    }
  }
  
  boost_poly_2
  to_boost_poly_2(const labeled_polygon_3& p) {
    boost_poly_2 pr;
    for (auto poly : p.vertices()) {
      boost::geometry::append(pr, boost::geometry::model::d2::point_xy<double>(poly.x, poly.y));
    }

    boost::geometry::interior_rings(pr).resize(p.holes().size());
    for (int i = 0; i < p.holes().size(); i++) {
      auto& h = p.holes()[i];
      for (auto p : h) {
        boost::geometry::append(pr, boost::geometry::model::d2::point_xy<double>(p.x, p.y), i);
      }
    }

    boost::geometry::correct(pr);
    
    return pr;
  }

  boost_poly_2 rotate_to_2D(const labeled_polygon_3& p) {
    const rotation r = rotate_from_to(p.normal(), point(0, 0, 1));
    return to_boost_poly_2(apply(r, p));
  }

  double area(const polygon_3& p) {
    auto td = rotate_to_2D(p);
    return bg::area(td);
  }

  bool
  add_to_existing_polys_as_hole(const polygon_3& next,
				std::vector<polygon_3>& polygons) {
    if (polygons.size() == 0) { return false; }

    const rotation r = rotate_from_to(next.normal(), point(0, 0, 1));
    boost_poly_2 next_boost_p = to_boost_poly_2(apply(r, next));

    int max_ind = -1;
    double outer_area_max = -1.0;
    cout << "outer max = " << outer_area_max << endl;
    for (unsigned i = 0; i < polygons.size(); i++) {
      const auto& current_poly = polygons[i];

      // cout << "Current polygon" << endl;
      // vtk_debug_polygon(current_poly);
      // cout << "Both polygons" << endl;
      // vtk_debug_polygons({current_poly, next});

      boost_poly_2 current_boost_p = to_boost_poly_2(apply(r, current_poly));
      double outer_area = area(build_clean_polygon_3(current_poly.vertices()));
      cout << "outer area = " << outer_area << endl;

      bool within_current = bg::within(next_boost_p, current_boost_p);
      cout << "within_current = " << within_current << endl;

      if (within_current && outer_area > outer_area_max) {
	outer_area_max = outer_area;
	max_ind = i;
	cout << "set max_ind = " << max_ind << endl;
      }

      cout << "outer max = " << outer_area_max << endl;
    }

    if (max_ind >= 0) {
      auto& current_poly = polygons[max_ind];
      unsigned old_num_holes = current_poly.holes().size();
      current_poly.add_hole(next.vertices());

      unsigned new_num_holes = polygons[max_ind].holes().size();

      DBG_ASSERT(new_num_holes == (old_num_holes + 1));

      cout << "ADDED TO POLYGON AS HOLE" << endl;

      return true;
    }

    return false;
  }

  std::vector<polygon_3>
  arrange_rings(const std::vector<std::vector<point>>& rings) {
    cout << "# of rings before arranging = " << rings.size() << endl;
    if (rings.size() == 0) { return {}; }

    vector<polygon_3> ring_polys;
    for (auto r : rings) {
      ring_polys.push_back(build_clean_polygon_3(r));
    }

    vector<polygon_3> polygons;
    while (ring_polys.size() > 0) {
      auto next_poly =
	max_element(begin(ring_polys), end(ring_polys),
		    [](const polygon_3& l, const polygon_3& r) {
		      return area(l) < area(r);
		    });

      DBG_ASSERT(next_poly != end(ring_polys));

      polygon_3 next_p = *next_poly;

      ring_polys.erase(next_poly);

      bool is_hole = add_to_existing_polys_as_hole(next_p, polygons);

      if (!is_hole) {
	polygons.push_back(next_p);
      }
    }

    cout << "# of polygons after arranging = " << polygons.size() << endl;

    return polygons;
  }
  
}
