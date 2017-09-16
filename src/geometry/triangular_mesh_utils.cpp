#include "geometry/triangular_mesh_utils.h"

namespace gca {

  oriented_polygon max_area_outline(const std::vector<index_t>& inds,
				    const triangular_mesh& m) {
    auto part_outlines = mesh_bounds(inds, m);

    DBG_ASSERT(part_outlines.size() > 0);

    oriented_polygon part_outline =
      *(max_element(begin(part_outlines), end(part_outlines),
		    [](const oriented_polygon& l,
		       const oriented_polygon& r)
      { return area(l) < area(r); }));

    return part_outline;
  }

  vector<oriented_polygon> mesh_bounds(const vector<index_t>& faces,
				       const triangular_mesh& mesh) {
    vector<oriented_polygon> ps;
    if (faces.size() == 0) {
      return ps;
    }
    point normal = mesh.face_orientation(faces.front());
    typedef pair<index_t, index_t> iline;
    vector<iline> tri_lines;
    for (auto i : faces) {
      auto t = mesh.triangle_vertices(i);
      tri_lines.push_back(iline(t.v[0], t.v[1]));
      tri_lines.push_back(iline(t.v[1], t.v[2]));
      tri_lines.push_back(iline(t.v[2], t.v[0]));
    }

    // TODO: Change to sort and count, maybe add to system/algorithm?
    vector<iline> no_ds;
    for (auto l : tri_lines) {
      int count = 0;
      for (auto r : tri_lines) {
	if ((l.first == r.first && l.second == r.second) ||
	    (l.first == r.second && l.second == r.first)) {
	  count++;
	}
      }
      DBG_ASSERT(count > 0);
      if (count == 1) {
	no_ds.push_back(l);
      }
    }

    vector<line> no_dups;
    for (auto l : no_ds) {
      no_dups.push_back(line(mesh.vertex(l.first), mesh.vertex(l.second)));
    }

    return unordered_segments_to_polygons(normal, no_dups);
  }

  vector<polygon_3> surface_boundary_polygons(const vector<index_t>& faces,
					      const triangular_mesh& mesh) {
    auto bounds = mesh_bounds(faces, mesh);

    vector<vector<point> > bound_rings;
    for (auto& bound : bounds) {
      bound_rings.push_back(bound.vertices());
    }

    std::vector<polygon_3> polys = arrange_rings(bound_rings);

    return polys;
  }

  polygon_3 surface_boundary_polygon(const vector<index_t>& faces,
				     const triangular_mesh& mesh) {
    auto ps = surface_boundary_polygons(faces, mesh);

    DBG_ASSERT(ps.size() == 1);

    return ps.front();
  }

  

  boost::optional<shared_edge>
  common_edge(const index_t l,
	      const index_t r,
	      const triangular_mesh& part) {
    auto tl = part.triangle_vertices(l);
    auto tr = part.triangle_vertices(r);
    std::vector<index_t> shared_verts;

    for (unsigned i = 0; i < 3; i++) {
      for (unsigned j = 0; j < 3; j++) {
	if (tl.v[i] == tr.v[j]) {
	  shared_verts.push_back(tl.v[i]);
	}
      }
    }

    DBG_ASSERT(shared_verts.size() < 3);

    if (shared_verts.size() == 2) {
      return shared_edge{l, r, edge(shared_verts[0], shared_verts[1])};
    }

    return boost::none;
  }
  
  std::vector<shared_edge> all_shared_edges(const std::vector<index_t>& l_faces,
					    const std::vector<index_t>& r_faces,
					    const triangular_mesh& part) {
    vector<shared_edge> edges;
    for (auto l : l_faces) {
      for (auto r : r_faces) {
	auto ce = common_edge(l, r, part);
	if (ce) {
	  edges.push_back(*ce);
	}
      }
    }
    return edges;
  }

  index_t non_edge_vertex_1(const shared_edge e, const triangular_mesh& m) {
    triangle_t t1 = m.triangle_vertices(e.triangle_1);
    vector<index_t> all_verts{t1.v[0], t1.v[1], t1.v[2]};
    vector<index_t> edge_verts{e.e.l, e.e.r};
    subtract(all_verts, edge_verts);

    DBG_ASSERT(all_verts.size() == 1);

    return all_verts.front();
  }


  index_t non_edge_vertex_2(const shared_edge e, const triangular_mesh& m) {
    triangle_t t1 = m.triangle_vertices(e.triangle_2);
    vector<index_t> all_verts{t1.v[0], t1.v[1], t1.v[2]};
    vector<index_t> edge_verts{e.e.l, e.e.r};
    subtract(all_verts, edge_verts);

    DBG_ASSERT(all_verts.size() == 1);

    return all_verts.front();
  }
  
  bool is_valley_edge(const shared_edge e,
		      const triangular_mesh& m) {
    point na = m.face_orientation(e.triangle_1);
    point pa = m.vertex(non_edge_vertex_1(e, m));
    point pb = m.vertex(non_edge_vertex_2(e, m));

    return dot((pb - pa), na) > 0.0;
  }

  bool angle_eps(const shared_edge e,
		 const triangular_mesh& m,
		 const double angle,
		 const double tol) {
    point n1 = m.face_orientation(e.triangle_1);
    point n2 = m.face_orientation(e.triangle_2);
    return angle_eps(n1, n2, angle, tol);
  }

  polygon_3 box_bound(const triangular_mesh& mesh) {
    box bb = mesh.bounding_box();

    point p1(bb.x_min, bb.y_min, bb.z_min);
    point p2(bb.x_max, bb.y_min, bb.z_min);
    point p3(bb.x_max, bb.y_max, bb.z_min);
    point p4(bb.x_min, bb.y_max, bb.z_min);

    vector<point> pts{p1, p2, p3, p4};
    polygon_3 poly = build_clean_polygon_3(pts);
    poly.correct_winding_order(point(0, 0, 1));

    return poly;
  }

  double angle_between_normals(const shared_edge e,
			       const triangular_mesh& m) {
    point n1 = m.face_orientation(e.triangle_1);
    point n2 = m.face_orientation(e.triangle_2);
    return angle_between(n1, n2);
  }

  // Move to triangular mesh utils
  std::vector<index_t>
  vertex_inds_on_surface(const std::vector<index_t>& s,
		      const triangular_mesh& m) {
    std::vector<index_t> pts;
    std::unordered_set<index_t> already_added;

    for (auto j : s) {
      auto t = m.triangle_vertices(j);
      for (unsigned i = 0; i < 3; i++) {
	index_t v = t.v[i];
	if (already_added.find(v) == end(already_added)) {
	  already_added.insert(v);
	  pts.push_back(v);
	}
      }
    }
    
    return pts;
  }
  
}
