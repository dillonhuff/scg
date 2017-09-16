#include "geometry/offset.h"
#include "geometry/polygon_3.h"
#include "geometry/rotation.h"
#include "geometry/vtk_debug.h"

namespace gca {

  polygon_3
  build_clean_polygon_3(const std::vector<point>& vertices,
			const std::vector<std::vector<point>>& hole_verts) {
    vector<point> outer_ring = clean_vertices(vertices);
    delete_antennas(outer_ring);

    // There is an occasional test failure here in simple box
    if (!(outer_ring.size() >= 3)) {
      cout << "ERROR: Outer ring size = " << outer_ring.size() << endl;
      vtk_debug_ring(outer_ring);

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


  // NOTE: These should be references
  polygon_3::polygon_3(const std::vector<point>& vertices,
		       const std::vector<std::vector<point>>& hole_verts,
		       const bool dummy_param) :
    outer_ring(vertices),
    inner_rings{} {

    for (auto h : hole_verts) {
      inner_rings.push_back(h);
    }
  }
  
  void check_simplicity(const polygon_3& p) {
    check_simplicity(static_cast<const std::vector<point>&>(p.vertices()));

    for (auto h : p.holes()) {
      check_simplicity(static_cast<const std::vector<point>&>(h));
    }
  }

  std::vector<point> project_points(const plane pl,
				    const std::vector<point>& pts) {
    vector<point> res_pts;
    for (auto p : pts) {
      res_pts.push_back(project(pl, p));
    }

    return res_pts;
  }
  
  polygon_3 project_onto(const plane p,
			 const polygon_3& poly) {

    vector<point> proj_outer = project(p, poly.vertices());

    vector<vector<point>> proj_holes;
    for (auto h : poly.holes()) {
      proj_holes.push_back(project(p, h));
    }

    polygon_3 l = build_clean_polygon_3(proj_outer, proj_holes);

    if (!(within_eps(angle_between(l.normal(), p.normal()), 0.0, 0.1))) {
      l.correct_winding_order(p.normal());
    }
    
    DBG_ASSERT(within_eps(angle_between(l.normal(), p.normal()), 0.0, 0.1));
    
    return l;
  }

  boost_poly_2
  to_boost_poly_2(const labeled_polygon_3& p) {
    boost_poly_2 pr;
    for (auto p : p.vertices()) {
      boost::geometry::append(pr, boost::geometry::model::d2::point_xy<double>(p.x, p.y));
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

  polygon_3
  to_polygon_3(const double z, const boost_poly_2& p) {
    vector<point> vertices;
    for (auto p2d : boost::geometry::exterior_ring(p)) {
      point pt(p2d.get<0>(), p2d.get<1>(), z);
      vertices.push_back(pt);
    }

    vector<vector<point>> holes;
    for (auto ir : boost::geometry::interior_rings(p)) {
      vector<point> hole_verts;
      for (auto p2d : ir) {
	point pt(p2d.get<0>(), p2d.get<1>(), z);
	hole_verts.push_back(pt);
      }
      holes.push_back(clean_vertices(hole_verts));
    }
    return build_clean_polygon_3(clean_vertices(vertices), holes);

  }

  boost_poly_2 rotate_to_2D(const labeled_polygon_3& p) {
    const rotation r = rotate_from_to(p.normal(), point(0, 0, 1));
    return to_boost_poly_2(apply(r, p));
  }

  double area(const polygon_3& p) {
    auto td = rotate_to_2D(p);
    return bg::area(td);
  }

  boost_multipoly_2
  planar_union_boost(const std::vector<polygon_3>& polys) {
    if (polys.size() == 0) { return {}; }

    //vtk_debug_polygons(polys);

    double level_z =
      max_distance_along(polys.front().vertices(), polys.front().normal());
    point n = polys.front().normal();

    cout << "n = " << n << endl;
    
    const rotation r = rotate_from_to(n, point(0, 0, 1));
    const rotation r_inv = inverse(r);

    cout << "# polys to union = " << polys.size() << endl;

    boost_multipoly_2 result;
    result.push_back(to_boost_poly_2(apply(r, polys.front())));

    for (unsigned i = 1; i < polys.size(); i++) {
      auto& s = polys[i];

      auto bp = to_boost_poly_2(apply(r, s));
      boost_multipoly_2 r_tmp = result;
      boost::geometry::clear(result);
      boost::geometry::union_(r_tmp, bp, result);
    }

    cout << "# polys in result = " << result.size() << endl;

    return result;
  }

  std::vector<polygon_3>
  planar_polygon_union(const point n,
		       const std::vector<polygon_3>& polys) {
    if (polys.size() == 0) { return {}; }

    //vtk_debug_polygons(polys);

    double level_z =
      max_distance_along(polys.front().vertices(), n); //polys.front().normal());

    cout << "In planar union" << endl;
    cout << "n = " << n << endl;
    // for (auto p : polys) {
    //   cout << "normal = " << p.normal() << endl;
    // }
    
    const rotation r = rotate_from_to(n, point(0, 0, 1));
    const rotation r_inv = inverse(r);

    cout << "# polys to union = " << polys.size() << endl;

    boost_multipoly_2 result;
    result.push_back(to_boost_poly_2(apply_no_check(r, polys.front())));

    for (unsigned i = 1; i < polys.size(); i++) {
      auto& s = polys[i];

      auto bp = to_boost_poly_2(apply_no_check(r, s));
      boost_multipoly_2 r_tmp = result;
      boost::geometry::clear(result);
      boost::geometry::union_(r_tmp, bp, result);
    }

    cout << "# polys in result = " << result.size() << endl;
    

    std::vector<polygon_3> res;
    for (auto& r : result) {
      boost::optional<polygon_3> lp =
	to_labeled_polygon_3_maybe(r_inv, level_z, r);

      if (lp) {
	check_simplicity(*lp);

	(*lp).correct_winding_order(n);
	res.push_back(*lp);
      }
    }

    cout << "# of polygon_3 in final output = " << res.size() << endl;

    return res;
    
  }
  
  // TODO: Use planar_union_boost as starting point here?
  std::vector<polygon_3>
  planar_polygon_union(const std::vector<polygon_3>& polys) {
    if (polys.size() == 0) { return {}; }
    point n = polys.front().normal();

    return planar_polygon_union(n, polys);
  }

  labeled_polygon_3
  convex_hull_2D(const std::vector<point>& pts,
		 const point n,
		 const double z_level) {
    const rotation r = rotate_from_to(n, point(0, 0, 1));
    const rotation r_inv = inverse(r);
    
    auto rotated_pts = apply(r, pts);

    boost_multipoint_2 mp;
    for (auto p : rotated_pts) {
      boost::geometry::append(mp, boost::geometry::model::d2::point_xy<double>(p.x, p.y));
    }

    boost_multipoint_2 res;
    boost::geometry::convex_hull(mp, res);

    vector<point> res_pts;
    for (auto p : res) {
      point pz(p.x(), p.y(), z_level);
      res_pts.push_back(times_3(r_inv, pz));
    }

    return build_clean_polygon_3(clean_vertices(res_pts));
  }

  box bounding_box(const polygon_3& p) {
    return bound_positions(static_cast<std::vector<point>>(p.vertices()));
  }

  polygon_3 project(const polygon_3& p, double z) {
    vector<point> pts;
    for (auto pt : p.vertices()) {
      pts.push_back(point(pt.x, pt.y, z));
    }

    vector<vector<point>> holes;
    for (auto h : p.holes()) {

      vector<point> hole;
      for (auto pt : h) {
	hole.push_back(point(pt.x, pt.y, z));
      }

      holes.push_back(hole);
    }

    return build_clean_polygon_3(pts, holes);
  }

  double area(const std::vector<polygon_3>& p) {
    double total_area = 0.0;

    for (auto poly : p) {
      double poly_area = area(poly);

      total_area += poly_area;
    }

    return total_area;
  }

  boost_multipoly_2
  to_boost_multipoly_2(const std::vector<polygon_3>& lines) {
    boost_multipoly_2 res;
    for (auto& pl : lines) {
      res.push_back(to_boost_poly_2(pl));
    }
    return res;
  }

  polygon_3 shift(const point p, const polygon_3& poly) {
    auto dr = shift(p, poly.vertices());

    vector<vector<point>> dh;
    for (auto h : poly.holes()) {
      dh.push_back(shift(p, h));
    }

    polygon_3 shifted = build_clean_polygon_3(dr, dh);
    return shifted;
  }

  std::vector<polygon_3> shift(const point p,
			       const std::vector<polygon_3>& polys) {
    vector<polygon_3> shifted;
    for (auto& poly : polys) {
      shifted.push_back(shift(p, poly));
    }
    return shifted;
  }

  polyline to_polyline(const polygon_3& p) {
    DBG_ASSERT(p.holes().size() == 0);

    auto v = p.vertices();
    v.push_back(p.vertices().front());
    return polyline(v);
  }

  std::vector<polyline> to_polylines(const polygon_3& p) {

    vector<polyline> lines;
    auto v = p.vertices();
    v.push_back(p.vertices().front());

    lines.push_back(polyline(v));

    for (auto& h : p.holes()) {
      lines.push_back(polyline(h));
    }

    return lines;
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

  labeled_polygon_3 dilate(const labeled_polygon_3& p, const double tol) {
    cout << "P normal = " << p.normal() << endl;
    
    auto dr = exterior_offset(clean_vertices(p.vertices()), tol);

    vector<vector<point>> dh;
    for (auto h : p.holes()) {
      auto offs = interior_offsets(clean_vertices(h), tol);
      for (auto off : offs) {
	dh.push_back(clean_vertices(off));
      }
    }

    labeled_polygon_3 poly =
      build_clean_polygon_3(clean_vertices(dr), dh);
    return poly;
  }

  polygon_3 shrink(const polygon_3& p, const double tol) {
    auto dr = interior_offset(p.vertices(), tol);

    vector<vector<point>> dh;
    for (auto h : p.holes()) {
      dh.push_back(exterior_offset(h, tol));
    }

    polygon_3 poly =
      build_clean_polygon_3(dr, dh);
    return poly;
  }

  boost::optional<labeled_polygon_3>
  shrink_optional(const labeled_polygon_3& p,
		  const double tol) {

    double level_z =
      max_distance_along(p.vertices(), p.normal());
    
    auto drs = interior_offsets(p.vertices(), tol);

    if (drs.size() == 0) { return boost::none; }

    if (drs.size() != 1) {
      //      vtk_debug_polygon(p);

      cout << "# of interior offsets = " << drs.size() << endl;

      // for (auto r : drs) {
      // 	vtk_debug_ring(r);
      // }

      //      DBG_ASSERT(false);
      return boost::none;
    }

    auto dr = drs.front();
    auto dr_pts = clean_vertices(dr);

    vector<polygon_3> new_holes;
    for (auto h : p.holes()) {
      auto h_clean = clean_vertices(exterior_offset(h, tol));
      if (h_clean.size() >= 3) {
	new_holes.push_back(build_clean_polygon_3(h_clean));
      }
    }

    //    vector<polygon_3> result_holes = planar_polygon_union(new_holes);

    const rotation r = rotate_from_to(p.normal(), point(0, 0, 1));
    const rotation r_inv = inverse(r);

    boost_poly_2 new_exterior_p(to_boost_poly_2(apply(r, build_clean_polygon_3(dr_pts))));

    boost_multipoly_2 union_of_holes = planar_union_boost(new_holes);

    boost_multipoly_2 shrunk_res;
    bg::difference(new_exterior_p, union_of_holes, shrunk_res);

    std::vector<polygon_3> res;
    for (auto r : shrunk_res) {
      polygon_3 lp = to_labeled_polygon_3(r_inv, level_z, r);

      check_simplicity(lp);

      lp.correct_winding_order(p.normal());
      res.push_back(lp);
    }

    if (res.size() != 1) { return boost::none; }

    return res.front();
    
    
    // if (bg::within(new_exterior_p, union_of_holes)) {
    //   return boost::none;
    // }

    // // vector<vector<point>> dh;
    // // for (auto h : result_holes) {
    // //   dh.push_back(h.vertices());
    // // }
    
    // if (dr_pts.size() >= 3) {
    //   labeled_polygon_3 poly(dr_pts, dh);
    //   return poly;
    // }

    // return boost::none;
  }

  labeled_polygon_3 smooth_buffer(const labeled_polygon_3& p,
				  const double tol) {
    return shrink(dilate(p, tol), tol);
  }
  
  std::vector<labeled_polygon_3>
  dilate_polygons(const std::vector<labeled_polygon_3>& polys, const double tol) {

    std::vector<labeled_polygon_3> dilated_polys;
    for (auto p : polys) {
      dilated_polys.push_back(dilate(p, tol));
    }

    return dilated_polys;
  }


  int curve_count(const polygon_3& f) {
    int count = 0;
    count += curve_count(f.vertices());

    for (auto& h : f.holes()) {
      int ch = curve_count(h);
      count += ch;

      cout << "Hole curve count = " << ch << endl;
      //vtk_debug_ring(h);
    }

    return count;
  }


  vector<polygon_3> from_boost_multipoly_2(const boost_multipoly_2& p,
					   const rotation& r,
					   const double z_level) {
    vector<polygon_3> polys;
    for (auto pl : p) {
      polys.push_back(apply(r, to_polygon_3(z_level, pl)));
    }
    return polys;
  }
  
  
  std::vector<polygon_3>
  polygon_intersection(const std::vector<polygon_3>& as,
		       const std::vector<polygon_3>& bs) {
    if (as.size() == 0 || bs.size() == 0) { return {}; }

    const rotation r = rotate_from_to(as.front().normal(), point(0, 0, 1));
    const rotation r_inv = inverse(r);
    double z_level = as.front().vertices().front().z;

    boost_multipoly_2 ap = to_boost_multipoly_2(r, as);
    boost_multipoly_2 bp = to_boost_multipoly_2(r, bs);
    boost_multipoly_2 cp;
    bg::intersection(ap, bp, cp);

    return from_boost_multipoly_2(cp, r, z_level);
  }

  std::vector<polygon_3>
  polygon_union(const std::vector<polygon_3>& as,
		const std::vector<polygon_3>& bs) {
    if (as.size() == 0 || bs.size() == 0) { return {}; }

    const rotation r = rotate_from_to(as.front().normal(), point(0, 0, 1));
    const rotation r_inv = inverse(r);
    double z_level = as.front().vertices().front().z;

    boost_multipoly_2 ap = to_boost_multipoly_2(r, as);
    boost_multipoly_2 bp = to_boost_multipoly_2(r, bs);
    boost_multipoly_2 cp;
    bg::union_(ap, bp, cp);

    return from_boost_multipoly_2(cp, r, z_level);
  }

  bool
  contains(const std::vector<polygon_3>& as,
	   const std::vector<polygon_3>& bs) {
    if (as.size() == 0 || bs.size() == 0) { return {}; }

    const rotation r = rotate_from_to(as.front().normal(), point(0, 0, 1));
    const rotation r_inv = inverse(r);
    double z_level = as.front().vertices().front().z;

    boost_multipoly_2 ap = to_boost_multipoly_2(r, as);
    boost_multipoly_2 bp = to_boost_multipoly_2(r, bs);

    return bg::within(bp, ap);
  }
  
  std::vector<polygon_3>
  polygon_difference(const std::vector<polygon_3>& as,
		     const std::vector<polygon_3>& bs) {
    if (as.size() == 0 || bs.size() == 0) { return {}; }

    const rotation r = rotate_from_to(as.front().normal(), point(0, 0, 1));
    const rotation r_inv = inverse(r);
    double z_level = as.front().vertices().front().z;

    boost_multipoly_2 ap = to_boost_multipoly_2(r, as);
    boost_multipoly_2 bp = to_boost_multipoly_2(r, bs);
    boost_multipoly_2 cp;
    bg::difference(ap, bp, cp);

    return from_boost_multipoly_2(cp, r, z_level);
  }
  
}
