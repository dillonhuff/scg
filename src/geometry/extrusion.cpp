#include "geometry/vtk_debug.h"
#include <boost/optional.hpp>

#include "geometry/extern_triangulate.h"
#include "geometry/extrusion.h"
#include "geometry/offset.h"
#include "geometry/rotation.h"
#include "geometry/triangular_mesh_utils.h"
#include "utils/relation.h"

#define ANSI_DECLARATORS

#define REAL double

#include <string>
#include <stdio.h>
#include <stdlib.h>

extern "C" {
#include "triangle_lib/triangle.h"
}

namespace gca {

  oriented_polygon
  oriented_polygon_for_index_polyline(const std::vector<point>& pts,
				      const index_poly& p,
				      const point n) {
    vector<point> poly_pts;
    for (unsigned in = 0; in < p.size(); in++) {
      if (in < p.size() - 1) {
	poly_pts.push_back(pts[p[in]]);
      } else if (p[in] != p[0]) {
	poly_pts.push_back(pts[p[in]]);
      }
    }
    return oriented_polygon(n, poly_pts);
  }
  
  oriented_polygon
  convert_index_poly(const std::vector<point>& pts,
		     const index_poly& p,
		     const double d,
		     const point e) {
    vector<point> poly_pts;
    for (auto i : p) {
      poly_pts.push_back(pts[i] + d*e.normalize());
    }
    return oriented_polygon(e, poly_pts);
  }

  triangular_mesh
  mesh_for_polys(const std::vector<oriented_polygon>& polys) {
    cout << "# of polygons = " << polys.size() << endl;

    std::vector<triangle> tris;
    for (auto p : polys) {
      concat(tris, vtk_triangulate_poly(p));
    }

    return make_mesh(tris, 0.01);
  }

  oriented_polygon
  side_rect(const std::vector<point>& pts,
	    const gca::edge e,
	    const double d0,
	    const double d1,
	    const point dir) {
    assert(d1 > d0);
    point i0 = d0*dir.normalize();
    point i1 = d1*dir.normalize();
    point top_left = pts[e.l] + i0;
    point top_right = pts[e.r] + i0;
    point bottom_left = pts[e.l] + i1;
    point bottom_right = pts[e.r] + i1;
    std::vector<point> poly_points{bottom_left, bottom_right, top_right, top_left};
    // TODO: Proper normal computation
    return oriented_polygon(dir, poly_points);
  }

  std::vector<gca::edge>
  index_poly_edges(const index_poly& p) {
    std::vector<gca::edge> edges;
    for (index_t i = 0; i < p.size(); i++) {
      index_t j = (i + 1) % p.size();
      edges.push_back(edge(p[i], p[j]));
    }
    return edges;
  }
  
  std::vector<oriented_polygon>
  side_polys(const std::vector<point>& pts,
	     const index_poly& p,
	     const double d0,
	     const double d1,
	     const point e) {
    auto edges = index_poly_edges(p);
    std::vector<oriented_polygon> side_rects;
    for (auto ed : edges) {
      side_rects.push_back(side_rect(pts, ed, d0, d1, e));
    }
    return side_rects;
  }

  template<typename T>
  std::vector<T>
  symmetric_difference(const std::vector<T>& a,
		       const std::vector<T>& b) {
    auto ac = a;
    subtract(ac, b);

    auto bc = b;
    subtract(bc, a);

    concat(ac, bc);

    return ac;
  }

  // Templatize and merge with collect polygon in triangle.cpp?
  index_poly collect_polygon(vector<gca::edge>& lines) {
    assert(lines.size() > 0);
    vector<index_t> points;
    vector<gca::edge> to_remove;
    points.push_back(lines.front().l);
    points.push_back(lines.front().r);
    lines.erase(lines.begin());
    unsigned i = 0;
    while (lines.size() > 0 && i < lines.size()) {
      if (lines[i].l == points.back()) {
    	if (lines[i].r == points.front()) {
    	  lines.erase(lines.begin() + i);	  
    	  return points;
    	}
    	points.push_back(lines[i].r);
    	lines.erase(lines.begin() + i);
    	i = 0;
      } else if (lines[i].r == points.back()) {
    	if (lines[i].l == points.front()) {
    	  lines.erase(lines.begin() + i);
    	  return points;
    	}
    	points.push_back(lines[i].l);
    	lines.erase(lines.begin() + i);
    	i = 0;
      } else {
    	i++;
      }
    }
    return points;
  }

  std::vector<index_poly>
  unordered_segments_to_index_polygons(std::vector<gca::edge>& lines) {
    vector<index_poly> ps;
    auto ls = lines;
    while (ls.size() > 0) {
      index_poly vertices = collect_polygon(ls);
      ps.push_back(vertices);
    }
    return ps;
  }

  relation<index_poly, index_t>
  build_endpoint_relation(const std::vector<index_poly>& plines,
			  const std::vector<index_t>& endpoints) {
    relation<index_poly, index_t> rel(plines, endpoints);
    for (auto i : rel.left_inds()) {
      const index_poly& p = rel.left_elem(i);
      for (auto j : rel.right_inds()) {
	index_t pt = rel.right_elem(j);
	if (p.front() == pt) {
	  rel.insert(i, j);
	}
	if (p.back() == pt) {
	  rel.insert(i, j);
	}
      }
    }
    return rel;
  }

  boost::optional<index_poly>
  merge_center(const index_poly& l, const index_poly& r) {
    if (l.back() == r.front()) {
      index_poly rest(begin(r) + 1, end(r));
      index_poly lc = l;
      concat(lc, rest);
      return lc;
    }
    return boost::none;
  }
  
  boost::optional<index_poly>
  merge_adjacent(const index_poly& l, const index_poly& r) {
    auto res = merge_center(l, r);
    if (res) { return *res; }
    res = merge_center(r, l);
    if (res) { return *res; }
    index_poly lc = l;
    reverse(begin(lc), end(lc));
    res = merge_center(lc, r);
    if (res) { return *res; }
    index_poly rc = r;
    reverse(begin(rc), end(rc));
    res = merge_center(l, rc);
    if (res) { return *res; }
    return boost::none;
  }

  bool
  try_to_merge_lines(std::vector<index_poly>& plines) {
    for (unsigned i = 0; i < plines.size(); i++) {
      index_poly* pi = &(plines[i]);
      for (unsigned j = 0; j < plines.size(); j++) {
	index_poly* pj = &(plines[j]);
	if (i != j) {
	  boost::optional<index_poly> res = merge_adjacent(*pi, *pj);
	  if (res) {
	    *pi = *res;
	    plines.erase(begin(plines) + j);
	    return true;
	  }
	}
      }
    }
    return false;
  }
  
  std::vector<index_poly>
  unordered_segments_to_index_polylines(std::vector<gca::edge>& lines) {
    assert(lines.size() > 0);
    vector<index_poly> plines;
    for (auto e : lines) {
      plines.push_back({e.l, e.r});
    }
    bool merged_one = true;
    while (merged_one) {
      merged_one = try_to_merge_lines(plines);
    }
    return plines;
  }
  
  triangular_mesh
  extrude_layers(const std::vector<point>& pts,
		 const std::vector<index_poly>& poly_layers,
		 const std::vector<double>& layer_depths,
		 const point extrude_dir) {
    assert(pts.size() > 0);
    assert(layer_depths.size() == poly_layers.size());
    std::vector<oriented_polygon> polys;

    double last_depth_offset = 0.0;
    std::vector<gca::edge> last_edges{};
    for (unsigned i = 0; i < layer_depths.size(); i++) {
      auto current_edges = index_poly_edges(poly_layers[i]);
      // Front polys
      auto front_edges = symmetric_difference(current_edges, last_edges);
      auto front_polys =
	unordered_segments_to_index_polygons(front_edges);

      for (auto p : front_polys) {
	polys.push_back(convert_index_poly(pts, p, last_depth_offset, extrude_dir));
      }


      // Sides
      concat(polys, side_polys(pts, poly_layers[i], last_depth_offset, last_depth_offset + layer_depths[i], extrude_dir));

      // Back
      if (i + 1 == layer_depths.size()) {
	// auto last_edges = index_poly_edges(poly_layers[i]);
	// auto last_polys = unordered_segments_to_index_polygons(last_edges);
	// assert(last_polys.size() == 1);
	oriented_polygon back =
	  convert_index_poly(pts, poly_layers[i], last_depth_offset + layer_depths[i], extrude_dir);
	vector<point> verts = back.vertices();
	reverse(begin(verts), end(verts));
	oriented_polygon back_rev(-1*extrude_dir, verts);
	polys.push_back(back_rev);
      }

      last_edges = current_edges;
      last_depth_offset += layer_depths[i];
    }
    return mesh_for_polys(polys);
  }

  triangular_mesh
  extrude(const extrusion& ext) {
    return extrude_layers(ext.pts, ext.poly_layers, ext.layer_depths, ext.extrude_dir);
  }


  // TODO: Clarify open vs. closed
  polyline
  to_polyline(const index_poly& poly,
	      const std::vector<point>& pts) {
    vector<point> res;
    for (auto i : poly) {
      res.push_back(pts[i]);
    }
    return res;
  }

  // TODO: Make this less hacky
  index_poly
  min_index_poly(const std::vector<point>& pts,
		 const std::vector<index_poly>& polys) {
    DBG_ASSERT(polys.size() > 0);
    return *(min_element(begin(polys), end(polys),
			 [pts](const index_poly& l, const index_poly& r)
			 { return pts[l.front()].z < pts[r.front()].z; }));
  }

  std::vector<gca::edge>
  index_poly_to_edges(const index_poly& p) {
    DBG_ASSERT(p.size() > 0);
    vector<gca::edge> edges;
    for (unsigned i = 0; i < p.size(); i++) {
      gca::edge e(p[i], p[(i + 1) % p.size()]);
      edges.push_back(e);
    }

    DBG_ASSERT(p.size() == edges.size());
    
    return edges;
  }

  std::vector<triangle> triangulate(const polygon_3& p) {
    const rotation r = rotate_from_to(p.normal(), point(0, 0, 1));
    const rotation r_inv = inverse(r);

    polygon_3 r_p = apply(r, p);

    DBG_ASSERT(angle_eps(r_p.normal(), point(0, 0, 1), 0.0, 1.0));
    
    std::vector<triangle> flat_3d_tris = triangulate_flat_3d(r_p);

    std::vector<triangle> res;
    for (auto t : flat_3d_tris) {
      res.push_back(apply(r_inv, t));
    }

    return res;
  }

  std::vector<triangle> build_top(const std::vector<triangle>& base_tris,
				  const point v) {
    std::vector<triangle> top;
    for (auto t : base_tris) {
      triangle top_t = triangle(-1*t.normal,
				v + t.v2,
				v + t.v1,
				v + t.v3);
      top.push_back(top_t);
    }

    return top;
  }

  std::vector<triangle> ring_sides(const std::vector<point>& poly, const point v) {
    vector<triangle> side_tris;
    
    for (unsigned i = 0; i < poly.size(); i++) {
      point p = poly[i];
      point q = poly[(i + 1) % poly.size()];

      point r = p + v;
      point s = q + v;

      vector<point> side_triangle_1{p, q, r};
      point n1 = ring_normal(side_triangle_1);
      triangle t1(n1,
		  side_triangle_1[0],
		  side_triangle_1[1],
		  side_triangle_1[2]);

      vector<point> side_triangle_2{r, q, s};
      point n2 = ring_normal(side_triangle_2);
      triangle t2(n2,
		  side_triangle_2[0],
		  side_triangle_2[1],
		  side_triangle_2[2]);
      
      //      polygon_3 side_poly({p, q, s, r});
      vector<triangle> tris{t1, t2}; // = triangulate(side_poly);
      concat(side_tris, tris);
    }

    return side_tris;
  }

  std::vector<triangle> build_sides(const polygon_3& poly, const point v) {
    vector<triangle> side_tris;
    concat(side_tris, ring_sides(poly.vertices(), v));

    for (auto h : poly.holes()) {
      concat(side_tris, ring_sides(h, v));
    }

    return side_tris;
  }

  triangular_mesh extrude(const polygon_3& p, const point v) {
    std::vector<triangle> base = triangulate(p);

    std::vector<triangle> top = build_top(base, v);

    concat(base, top);

    std::vector<triangle> sides = build_sides(p, v);
    cout << "# of triangles in sides = " << sides.size() << endl;
    concat(base, sides);

    return make_mesh(base, 0.0001);
  }

  polygon_3 base_polygon(const std::vector<index_t>& surf,
			 const triangular_mesh& part) {
    // auto bounds = mesh_bounds(surf, part);

    // vector<vector<point> > bound_rings;
    // for (auto& bound : bounds) {
    //   bound_rings.push_back(bound.vertices());
    // }

    vector<polygon_3> polys =
      surface_boundary_polygons(surf, part); //arrange_rings(bound_rings);

    if (polys.size() != 1) {
      cout << "ERROR: polys.size() = " << polys.size() << endl;
      vtk_debug_highlight_inds(surf, part);
      for (auto& p : polys) {
	vtk_debug_polygon(p);
      }
      vtk_debug_polygons(polys);

      DBG_ASSERT(polys.size() == 1);
    }

    return polys.front();
  }

  std::vector<triangle>
  flipped_surface_triangles(const std::vector<index_t>& surf,
			    const triangular_mesh& part) {
    vector<triangle> tris;
    for (auto i : surf) {
      triangle_t tvs = part.triangle_vertices(i);

      triangle flipped_tri(-1*part.face_orientation(i),
			   part.vertex(tvs.v[1]),
			   part.vertex(tvs.v[0]),
			   part.vertex(tvs.v[2]));

      tris.push_back(flipped_tri);
    }

    return tris;
  }
  
  triangular_mesh
  extrude_surface_negative(const std::vector<index_t>& surf,
			   const triangular_mesh& part,
			   const point v,
			   const double length) {

    polygon_3 base_poly = base_polygon(surf, part);

    std::vector<triangle> base = flipped_surface_triangles(surf, part);

    std::vector<triangle> top = build_top(base, length*v);

    concat(base, top);

    std::vector<triangle> sides = build_sides(base_poly, length*v);
    cout << "# of triangles in sides = " << sides.size() << endl;
    concat(base, sides);

    //vtk_debug_triangles(base);

    return make_mesh(base, 0.0001).flip_winding_order();
  }


}
