#include <cassert>

#include "geometry/offset.h"
#include "geometry/ring.h"
#include "geometry/rotation.h"

#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>

#include "geometry/vtk_debug.h"

namespace gca {

  typedef boost::shared_ptr<Polygon_2> PolygonPtr ;
  typedef std::vector<PolygonPtr> PolygonPtrVector ;
  typedef CGAL::Straight_skeleton_2<K> Ss ;

  typedef CGAL::Polygon_with_holes_2<K> PolygonWithHoles ;
  typedef boost::shared_ptr<PolygonWithHoles> PolygonWithHolesPtr ;
  typedef std::vector<PolygonWithHolesPtr> PolygonWithHolesPtrVector;
  
  typedef boost::shared_ptr<Ss> SsPtr ;  

  void set_orientation(Polygon_2& to_offset) {
    if (!(to_offset.is_simple())) {
      DBG_ASSERT(false);
    }
      
    
    if (to_offset.orientation() == CGAL::CLOCKWISE) {
      to_offset.reverse_orientation();
    }

    DBG_ASSERT(to_offset.orientation() == CGAL::COUNTERCLOCKWISE);
    
  }

  void set_orientation_clockwise(Polygon_2& to_offset) {
    if (!(to_offset.is_simple())) {
      DBG_ASSERT(false);
    }
      
    
    if (to_offset.orientation() == CGAL::COUNTERCLOCKWISE) {
      to_offset.reverse_orientation();
    }

    DBG_ASSERT(to_offset.orientation() == CGAL::CLOCKWISE);
    
  }
  
  Polygon_2
  CGAL_polygon_for_oriented_polygon(const oriented_polygon& p) {
    Polygon_2 out;
    for (auto p : p.vertices()) {
      out.push_back(Point(p.x, p.y));
    }

    return out;
  }

  Polygon_2
  CGAL_polygon_for_points(const std::vector<point>& pts) {
    Polygon_2 out;
    for (auto p : pts) {
      out.push_back(Point(p.x, p.y));
    }

    return out;
  }

  vector<point>
  ring_for_CGAL_polygon(const Polygon_2& off_p,
			const double z) {
    vector<point> res_pts;
    for (auto it = CGAL::CGAL_SS_i::vertices_begin(off_p);
	 it != CGAL::CGAL_SS_i::vertices_end(off_p); ++it) {
      Point vert = *it;
      res_pts.push_back(point(vert.x(), vert.y(), z));
    }

    return res_pts;
  }

  boost_poly_2 CGAL_polygon_for_boost_poly(const Polygon_2& poly) {
    double z_level = 0;
    vector<point> pts = ring_for_CGAL_polygon(poly, z_level);
    return to_boost_poly_2(build_clean_polygon_3(pts));
  }

  oriented_polygon
  oriented_polygon_for_CGAL_polygon(const Polygon_2& off_p,
				    const double z,
				    const point n) {
    vector<point> res_pts;
    for (auto it = CGAL::CGAL_SS_i::vertices_begin(off_p);
	 it != CGAL::CGAL_SS_i::vertices_end(off_p); ++it) {
      Point vert = *it;
      res_pts.push_back(point(vert.x(), vert.y(), z));
    }
    return oriented_polygon(n, res_pts);
  }

  std::vector<oriented_polygon> exterior_offset(const oriented_polygon& q,
						const double inc) {
    DBG_ASSERT(q.vertices().size() > 0);

    oriented_polygon p;
    if (signed_area(q) < 0) {
      p = q;
    } else {
      vector<point> pts = q.vertices();
      reverse(begin(pts), end(pts));
      p = oriented_polygon(q.normal(), pts);
    }

    DBG_ASSERT(p.vertices().size() > 0);
    
    double z_va = p.vertices().front().z;
    Polygon_2 out;
    for (auto p : p.vertices()) {
      out.push_back(Point(p.x, p.y));
    }

    if (!(out.is_simple())) {
      vtk_debug_polygon(p);
      DBG_ASSERT(false);
    }
    
    if (out.orientation() == CGAL::CLOCKWISE) {
      out.reverse_orientation();
    }

    DBG_ASSERT(out.orientation() == CGAL::COUNTERCLOCKWISE);

    DBG_ASSERT(out.is_simple());

    PolygonPtrVector inner_offset_polygons =
      CGAL::create_exterior_skeleton_and_offset_polygons_2(inc, out);
    
    vector<oriented_polygon> results;
    for (auto off_ptr : inner_offset_polygons) {
      Polygon_2 off_p = *off_ptr;
      auto op = oriented_polygon_for_CGAL_polygon(off_p, z_va, p.normal());
      results.push_back(op);
    }
    return results;
  }

  std::vector<oriented_polygon> interior_offset(const oriented_polygon& p,
						const double inc) {

    DBG_ASSERT(p.vertices().size() > 0);

    check_simplicity(p);

    double z_va = p.vertices().front().z;
    Polygon_2 out;
    for (auto p : p.vertices()) {
      out.push_back(Point(p.x, p.y));
    }

    if (!(out.is_simple())) {
      vtk_debug_polygon(p);
      DBG_ASSERT(false);
    }
    
    if (out.orientation() == CGAL::CLOCKWISE) {
      out.reverse_orientation();
    }

    DBG_ASSERT(out.orientation() == CGAL::COUNTERCLOCKWISE);

    //    cout << "Starting to compute offset skeleton" << endl;
    SsPtr ss = CGAL::create_interior_straight_skeleton_2(out);
    //    cout << "Done computing to offset skeleton" << endl;

    if (!ss) {
      cout << "straight skeleton is null" << endl;
      for (unsigned i = 0; i < p.vertices().size(); i++) {
	point pt = p.vertices()[i];
	unsigned i1 = (i + 1) % p.vertices().size();
	point pt1 = p.vertices()[i1];

	point dir = (pt1 - pt).normalize();
	
	cout << pt << "         with direction = " << dir << endl;
      }

      vtk_debug_ring(p.vertices());
      DBG_ASSERT(false);
    }

    //    cout << "Starting to compute the interior offset" << endl;
    PolygonPtrVector inner_offset_polygons =
      CGAL::create_offset_polygons_2<Polygon_2>(inc, *ss);
      //  CGAL::create_interior_skeleton_and_offset_polygons_2(inc, out);
    //    cout << "Done computing interior offset" << endl;

    vector<oriented_polygon> results;
    for (auto off_ptr : inner_offset_polygons) {
      Polygon_2 off_p = *off_ptr;
      auto op = oriented_polygon_for_CGAL_polygon(off_p, z_va, p.normal());
      results.push_back(op);
    }
    return results;
  }

  void check_simplicity(const oriented_polygon& p) {
    if (!CGAL_polygon_for_oriented_polygon(p).is_simple()) {
      vtk_debug_polygon(p);
      DBG_ASSERT(false);
    }
  }

  void check_simplicity(const std::vector<point>& rpts) {
    point norm = ring_normal(rpts);

    if (within_eps(norm.len(), 0.0, 0.0000001)) {
      cout << "Error: Degenerate ring!" << endl;
      cout << "# of points in ring = " << rpts.size() << endl;
      for (auto pt : rpts) {
	cout << pt << endl;
      }
      vtk_debug_ring(rpts);
      DBG_ASSERT(false);
    }

    const gca::rotation r = rotate_from_to(norm, point(0, 0, 1));
    auto pts = apply(r, rpts);

    if (!(no_duplicate_points(pts, 0.0001))) {
      vtk_debug_ring(pts);
      DBG_ASSERT(false);
    }
    
    if (!CGAL_polygon_for_points(pts).is_simple()) {
      cout << "Non simple ring!" << endl;
      for (auto p : pts) {
	cout << p.x << " , " << p.y << " , " << p.z << endl;
      }
      vtk_debug_ring(pts);
      DBG_ASSERT(false);
    }
  }

  vector<point> clean_ring_for_offsetting_no_fail(const std::vector<point>& ring) {
    vector<point> r_pts = clean_vertices_within_eps(ring, 0.005, 0.0000001);

    if (r_pts.size() < 3) { return r_pts; }

    delete_antennas(r_pts);

    return r_pts;
  }
  
  vector<point> clean_ring_for_offsetting(const std::vector<point>& ring) {
    vector<point> r_pts = clean_vertices_within_eps(ring, 0.005, 0.0000001);

    if (r_pts.size() < 3) { return r_pts; }

    auto old_rpts = r_pts;

    delete_antennas(r_pts);

    if (r_pts.size() < 3) {
      cout << "Before cleaning" << endl;
      cout << "# of points = " << ring.size() << endl;

      vtk_debug_ring(ring);

      cout << "After cleaning" << endl;
      cout << "# of points = " << old_rpts.size() << endl;

      vtk_debug_ring(old_rpts);

      cout << "# of points after deleting antennas = " << r_pts.size() << endl;

      DBG_ASSERT(false);
    }

    if (has_antenna(r_pts)) {
      DBG_ASSERT(false);
    }

    check_simplicity(r_pts);

    return r_pts;
  }

  polygon_3 clean_polygon_for_offsetting(const polygon_3& poly) {
    vector<point> cleaned_outer = clean_ring_for_offsetting(poly.vertices());
    vector<vector<point>> cleaned_holes;
    for (auto h : poly.holes()) {
      cleaned_holes.push_back(clean_ring_for_offsetting(h));
    }

    return build_clean_polygon_3(cleaned_outer, cleaned_holes);
  }

  boost::optional<polygon_3>
  clean_polygon_for_offsetting_maybe(const polygon_3& poly) {
    vector<point> cleaned_outer = clean_ring_for_offsetting_no_fail(poly.vertices());

    if (cleaned_outer.size() < 3) { return boost::none; }

    vector<vector<point>> cleaned_holes;
    for (auto h : poly.holes()) {
      auto cl_h = clean_ring_for_offsetting_no_fail(h);
      if (cl_h.size() >= 3) {
	cleaned_holes.push_back(cl_h);
      }
    }

    return build_clean_polygon_3(cleaned_outer, cleaned_holes);
  }
  
  std::vector<point> exterior_offset(const std::vector<point>& pts,
				     const double tol) {

    DBG_ASSERT(pts.size() > 2);
    
    point n(0, 0, 1);
    const rotation r = rotate_from_to(ring_normal(pts), n);
    const rotation r_inv = inverse(r);
    auto r_pts = apply(r, pts);

    auto original_rpts = r_pts;

    r_pts = clean_vertices_within_eps(r_pts, 0.005, 0.0000001);

    auto old_rpts = r_pts;

    delete_antennas(r_pts);

    if (r_pts.size() < 3) {
      cout << "Before cleaning" << endl;
      cout << "# of points = " << pts.size() << endl;

      vtk_debug_ring(pts);

      cout << "After cleaning" << endl;
      cout << "# of points = " << old_rpts.size() << endl;

      vtk_debug_ring(old_rpts);

      cout << "# of points after deleting antennas = " << r_pts.size() << endl;

      DBG_ASSERT(false);
    }

    if (has_antenna(r_pts)) {
      DBG_ASSERT(false);
    }

    check_simplicity(r_pts);

    
    auto res = exterior_offset(oriented_polygon(n, r_pts), tol);

    if (res.size() != 2) {
      cout << "Wrong number of exterior offsets!" << endl;
      cout << "# of exterior offsets = " << res.size() << endl;
      cout << "tol = " << tol << endl;

      cout << "# of points in ring = " << r_pts.size() << endl;
      cout << "vector<point> pts{" << endl;
      for (auto p : r_pts) {
	cout << "point(" << p.x << ", " << p.y << ", " << p.z << ")" << ", " << endl;
      }
      cout << "};" << endl;



      vtk_debug_ring(r_pts);
      for (auto r : res) {
	vtk_debug_ring(r.vertices());
      }

      cout << "res" << endl;
      for (auto res_p : res) {
	vtk_debug_polygon(res_p);
      }

      DBG_ASSERT(res.size() == 2);
    }

    auto rpoly = res[1];
    auto rpts = apply(r_inv, rpoly.vertices());

    correct_winding_order(rpts, ring_normal(pts));

    return rpts;
  }

  std::vector<point> interior_offset(const std::vector<point>& pts, const double tol) {
    point n(0, 0, 1);
    const rotation r = rotate_from_to(ring_normal(pts), n);
    const rotation r_inv = inverse(r);
    auto res = interior_offset(oriented_polygon(n, apply(r, pts)), tol);

    DBG_ASSERT(res.size() == 1);

    auto rpoly = res.front();
    auto rpts = apply(r_inv, rpoly.vertices());

    correct_winding_order(rpts, ring_normal(pts));

    return rpts;
  }

  std::vector<std::vector<point>>
  interior_offsets(const std::vector<point>& pts,
		   const double tol) {
    point n(0, 0, 1);
    const rotation r = rotate_from_to(ring_normal(pts), n);
    const rotation r_inv = inverse(r);

    cout << "# of points = " << pts.size() << endl;
    cout << "ring normal = " << ring_normal(pts) << endl;
    cout << "tol         = " << tol << endl;
    //vtk_debug_ring(pts);

    auto r_pts = apply(r, pts);

    check_simplicity(r_pts);

    //vtk_debug_ring(r_pts);
    
    auto res = interior_offset(oriented_polygon(n, r_pts), tol);


    vector<vector<point>> result_pts;
    for (auto rpoly : res) {
      auto rpts = apply(r_inv, rpoly.vertices());
      correct_winding_order(rpts, ring_normal(pts));
      result_pts.push_back(rpts);
    }

    return result_pts;
  }
  

  polygon_3 exterior_offset_flat(const polygon_3& poly,
				 const double d) {
    DBG_ASSERT(angle_eps(poly.normal(), point(0, 0, 1), 0.0, 1.0));
    DBG_ASSERT(poly.vertices().size() >= 3);

    double z_level = poly.vertices().front().z;

    Polygon_2 outer;
    for (auto p : poly.vertices()) {
      outer.push_back(Point(p.x, p.y));
    }

    set_orientation(outer);

    PolygonPtrVector offset_poly_with_holes =
      create_exterior_skeleton_and_offset_polygons_2(d, outer);

    DBG_ASSERT(offset_poly_with_holes.size() > 1);

    offset_poly_with_holes.erase(begin(offset_poly_with_holes));

    vector<vector<point>> pts;

    for (PolygonPtr& p : offset_poly_with_holes) {
      pts.push_back(clean_ring_for_offsetting(ring_for_CGAL_polygon(*p, z_level)));
    }

    vector<polygon_3> result_polys = arrange_rings(pts);

    DBG_ASSERT(result_polys.size() == 1);

    polygon_3 result_poly = result_polys.front();

    // // NOTE: Assumes the 2nd polygon returned by CGAL is always the exteriormost
    // // ring of the exterior offset after the containing rectangle
    // vector<point> outer_ring =
    //   ring_for_CGAL_polygon(*(offset_poly_with_holes.front()), z_level);

    boost_poly_2 outer_poly = to_boost_poly_2(result_poly);

    boost_multipoly_2 holes_to_subtract;
    for (auto h : poly.holes()) {
      polygon_3 hole_polygon = build_clean_polygon_3(h);
      vector<polygon_3> hole_interior_offsets =
	interior_offset({hole_polygon}, d);

      for (polygon_3& hole_poly : hole_interior_offsets) {
	DBG_ASSERT(hole_poly.holes().size() == 0);

	holes_to_subtract.push_back(to_boost_poly_2(hole_poly));
      }
    }

    // for (unsigned i = 1; i < offset_poly_with_holes.size(); i++) {
    //   vector<point> h_ring =
    // 	ring_for_CGAL_polygon(*(offset_poly_with_holes[i]), z_level);

    //   holes_to_subtract.push_back(to_boost_poly_2(h_ring));
    // }

    boost_multipoly_2 ext_offset;
    bg::difference(outer_poly, holes_to_subtract, ext_offset);

    if (!(ext_offset.size() == 1)) {
      cout << "OFFSETTING POLYGON BY " << d << endl;
      cout << "ext_offset.size() = " << ext_offset.size() << endl;
      vtk_debug_polygon(poly);

      cout << "OFFSET RESULTS" << endl;
      for (auto r : offset_poly_with_holes) {
	vector<point> outer_ring =
	  ring_for_CGAL_polygon(*r, z_level);
	vtk_debug_ring(outer_ring);

      }

      DBG_ASSERT(ext_offset.size() == 1);
    }

    polygon_3 exterior_poly = to_polygon_3(z_level, ext_offset.front());

    return exterior_poly;
    
    // TODO: Convert result of outer exterior offset into a boost polygon,
    // do interior offsets of holes, convert results to boost polygons and
    // then subtract them from original result

    //DBG_ASSERT(false);

    // vector<polygon_3> result_polys;

    // for (PolygonWithHolesPtr& p : offset_poly_with_holes) {
    //   PolygonWithHoles& pg = *p;
    //   Polygon_2 outer = pg.outer_boundary();
    //   vector<point> outer_pts = ring_for_CGAL_polygon(outer, z_level);

    //   vector<vector<point>> holes;
    //   for (auto it = pg.holes_begin(); it != pg.holes_end(); ++it) {
    // 	Polygon_2& cgal_hole = *it;
    // 	holes.push_back(ring_for_CGAL_polygon(cgal_hole, z_level));
    //   }

    //   result_polys.push_back(polygon_3(outer_pts, holes));
    // }

    // return result_polys;
  }

  polygon_3 exterior_offset(const polygon_3& poly,
  			    const double d) {
    point n(0, 0, 1);
    const rotation r = rotate_from_to(poly.normal(), n);
    const rotation r_inv = inverse(r);

    polygon_3 r_poly = apply(r, poly);

    check_simplicity(r_poly);

    polygon_3 res = exterior_offset_flat(r_poly, d);

    polygon_3 final_res = apply(r_inv, res);
    final_res.correct_winding_order(poly.normal());

    return final_res;
  }
  
  std::vector<polygon_3> exterior_offset(const std::vector<polygon_3>& polys,
					 const double d) {
    vector<polygon_3> exter_offsets;
    for (auto poly : polys) {
      exter_offsets.push_back(exterior_offset(poly, d));
    }

    return planar_polygon_union(exter_offsets);
  }

  std::vector<polygon_3> interior_offsets_flat(const polygon_3& ply,
					       const double d) {

    DBG_ASSERT(angle_eps(ply.normal(), point(0, 0, 1), 0.0, 1.0));
    DBG_ASSERT(ply.vertices().size() >= 3);

    polygon_3 poly = clean_polygon_for_offsetting(ply);

    double z_level = poly.vertices().front().z;

    Polygon_2 outer;
    for (auto p : poly.vertices()) {
      outer.push_back(Point(p.x, p.y));
    }

    set_orientation(outer);

    PolygonWithHoles to_offset( outer );

    for (auto h : poly.holes()) {
      Polygon_2 hole;

      for (auto pt : h) {
	hole.push_back( Point(pt.x, pt.y) );
      }

      set_orientation_clockwise(hole);
      to_offset.add_hole( hole );
    }

    PolygonPtrVector offset_poly_with_holes =
      CGAL::create_interior_skeleton_and_offset_polygons_2(d, to_offset);

    vector<vector<point>> pts;

    for (PolygonPtr& p : offset_poly_with_holes) {
      auto r_new =
	clean_ring_for_offsetting_no_fail(ring_for_CGAL_polygon(*p, z_level));

      if (r_new.size() >= 3) {
	check_simplicity(r_new);
	pts.push_back(r_new);
      }
    }

    vector<polygon_3> result_polys = arrange_rings(pts);

    return result_polys;
  }
  
  std::vector<polygon_3> interior_offsets(const polygon_3& poly,
					  const double d) {
    point n(0, 0, 1);
    const rotation r = rotate_from_to(poly.normal(), n);
    const rotation r_inv = inverse(r);

    polygon_3 r_poly = apply(r, poly);

    check_simplicity(r_poly);

    auto res = interior_offsets_flat(r_poly, d);

    // NOTE: Need to add z value adjustment?
    vector<polygon_3> result_pts;
    for (auto res_poly : res) {
      auto rpts = apply(r_inv, res_poly);
      rpts.correct_winding_order(poly.normal());

      result_pts.push_back(rpts);
    }

    return result_pts;
  }

  std::vector<polygon_3> interior_offset(const std::vector<polygon_3>& polys,
					 const double d) {
    vector<polygon_3> inter_offsets;
    for (auto poly : polys) {
      concat(inter_offsets, interior_offsets(poly, d));
    }
    return inter_offsets;
  }
  
}
