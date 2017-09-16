#include "geometry/surface.h"
#include "geometry/triangular_mesh_utils.h"
#include "geometry/vtk_debug.h"

namespace gca {

  bool surfaces_share_edge(const surface& l,
			   const surface& r) {
    auto ind1 = l.index_list();
    auto ind2 = r.index_list();
    return share_edge(ind1, ind2, l.get_parent_mesh());
  }
  
  bool surfaces_share_edge(const unsigned i,
			   const unsigned j,
			   const std::vector<surface>& surfaces) {
    auto ind1 = surfaces[i].index_list();
    auto ind2 = surfaces[j].index_list();
    return share_edge(ind1, ind2, surfaces[i].get_parent_mesh());
  }

  bool surfaces_share_edge(const unsigned i,
			   const unsigned j,
			   const std::vector<surface*>& surfaces) {
    auto ind1 = surfaces[i]->index_list();
    auto ind2 = surfaces[j]->index_list();
    return share_edge(ind1, ind2, surfaces[i]->get_parent_mesh());
  }
  
  std::vector<index_t> surface_vertexes(const surface& s) {
    vector<index_t> inds;
    for (auto i : s.index_list()) {
      triangle_t t = s.get_parent_mesh().triangle_vertices(i);
      inds.push_back(t.v[0]);
      inds.push_back(t.v[1]);
      inds.push_back(t.v[2]);
    }
    sort(begin(inds), end(inds));
    return inds;
  }

  bool orthogonal_flat_surfaces(const surface* l, const surface* r) {
    point l_orient = l->face_orientation(l->front());
    point r_orient = r->face_orientation(r->front());
    double theta = angle_between(l_orient, r_orient);
    return within_eps(theta, 90, 0.1);
  }

  bool parallel_flat_surfaces(const surface* l, const surface* r) {
    point l_orient = l->face_orientation(l->front());
    point r_orient = r->face_orientation(r->front());
    double theta = angle_between(l_orient, r_orient);
    return within_eps(theta, 180, 0.1);
  }

  std::vector<surface> outer_surfaces(const triangular_mesh& part) {
    //    cout << "# of triangles in the part = " << part.face_indexes().size() << endl;

    auto const_orient_face_indices = const_orientation_regions(part);
    vector<surface> surfaces;

    //    cout << "# of const orientation regions = " << const_orient_face_indices.size() << endl;
    for (auto f : const_orient_face_indices) {
      surface s(&part, f);

      // cout << "Region normal            = " << normal(s) << endl;
      // cout << "# of triangles in region = " << s.index_list().size() << endl;

      //vtk_debug_highlight_inds(f, part);

      DBG_ASSERT(f.size() > 0);

      bool is_outer = is_outer_surface(f, part);

      //      cout << "Is outer                 = " << is_outer << endl;

      if (is_outer) {
	surfaces.push_back(surface(&part, f));
      }
    }

    return surfaces;
  }

  surface merge_surfaces(const std::vector<surface>& surfaces) {
    DBG_ASSERT(surfaces.size() > 0);
    vector<index_t> inds;
    for (auto s : surfaces) {
      concat(inds, s.index_list());
    }
    return surface(&(surfaces.front().get_parent_mesh()), inds);
  }

  void remove_contained_surfaces(const std::vector<surface>& stable_surfaces,
				 std::vector<surface>& surfaces_to_cut) {
    vector<index_t> stable_surface_inds;
    for (auto s : stable_surfaces) {
      concat(stable_surface_inds, s.index_list());
    }
    sort(begin(stable_surface_inds), end(stable_surface_inds));

    delete_if(surfaces_to_cut,
	      [&stable_surface_inds](const surface& s)
	      { return s.contained_by_sorted(stable_surface_inds); });
  }

  std::vector<surface>
  merge_surface_groups(const std::vector<surface>& surfs,
		       const std::vector<std::vector<unsigned> >& groups) {
    std::vector<surface> sfs;
    for (auto group : groups) {
      vector<surface> sg = select_indexes(surfs, group);
      DBG_ASSERT(sg.size() == group.size());
      sfs.push_back(merge_surfaces(sg));
    }
    return sfs;
  }

  // TODO: Eventually identify vertical holes and compensate for them
  // TODO: Should this compensate for the fact that direction is not
  //
  bool vertical_contained_by(const surface& maybe_contained,
			     const surface& maybe_container) {
    // vector<oriented_polygon> maybe_contained_outlines =
    //   mesh_bounds(maybe_contained.index_list(), maybe_contained.get_parent_mesh());
    // if (maybe_contained_outlines.size() !=  2) {
    //   cout << "More than 2 outlines" << endl;
    //   return false;
    // }
    auto contained_outline =
      max_area_outline(maybe_contained.index_list(), maybe_contained.get_parent_mesh()); //maybe_contained_outlines.front();
    
    // vector<oriented_polygon> maybe_container_outlines =
    //   mesh_bounds(maybe_container.index_list(), maybe_container.get_parent_mesh());
    // if (maybe_container_outlines.size() != 2) {
    //   cout << "More than 2 outlines" << endl;
    //   return false;
    // }
    auto container_outline = max_area_outline(maybe_container.index_list(), maybe_container.get_parent_mesh()); //maybe_container_outlines.front();

    return contains(container_outline, contained_outline);
  }

  boost::optional<surface>
  part_outline_surface(std::vector<surface>* surfaces_to_cut,
		       const point n) {
    cout << "Computing part outline surface" << endl;
    vector<surface> vertical_surfs =
      select(*surfaces_to_cut,
    	     [n](const surface& s)
    	     { return s.orthogonal_to(n, 0.01); });

    cout << "# of vertical surfaces in " << n << " = " << vertical_surfs.size() << endl;

    vector<vector<unsigned>> merge_groups =
      connected_components_by(vertical_surfs, [](const surface& l, const surface& r)
    			      { return surfaces_share_edge(l, r); });

    if (merge_groups.size() == 1) {
      return merge_surfaces(vertical_surfs);
    } else {
      cout << "More than 1 vertical component" << endl;
      vector<surface> merged = merge_surface_groups(vertical_surfs, merge_groups);
      vector<surface> outer_surfs =
	partial_order_maxima(merged, [](const surface& l, const surface& r)
			     { return vertical_contained_by(l, r); });

      cout << "# of merged surfaces = " << merged.size() << endl;
      cout << "# of outer surfaces = " << outer_surfs.size() << endl;

      if (outer_surfs.size() == 1) {
	return outer_surfs.front();
      } else {
	return boost::none;
      }
    }

    return boost::none;
  }

  boost::optional<surface>
  part_outline_surface(const triangular_mesh& m,
		       const point n) {
    std::vector<surface> vertical_surfs =
      connected_vertical_surfaces(m, n);
    boost::optional<surface> outline =
      part_outline_surface(&vertical_surfs, n);
    return outline;
    // vector<surface> surfs = surfaces_to_cut(m);
    // return part_outline_surface(&surfs, n);
  }

  // TODO: Need to add normal vectors, how to match this with
  // the code in make_fixture_plan?
  boost::optional<oriented_polygon>
  part_outline(std::vector<surface>* surfaces_to_cut) {
    point n(0, 0, 1);

    auto m = part_outline_surface(surfaces_to_cut, n);
    if (m) {
      vector<oriented_polygon> outlines =
	mesh_bounds((*m).index_list(), (*m).get_parent_mesh());
      if (outlines.size() == 2) {
	vector<surface> vertical_surfs{*m};
	remove_contained_surfaces(vertical_surfs, *surfaces_to_cut);
	return outlines.front();
      }
    }

    return boost::none;
  }

  std::vector<surface> surfaces_to_cut(const triangular_mesh& part) {
    auto inds = part.face_indexes();

    return surfaces_to_cut(inds, part);
  }
  
  std::vector<surface> surfaces_to_cut(const std::vector<index_t>& indexes,
				       const triangular_mesh& part) {
    vector<index_t> inds = indexes;
    double normal_degrees_delta = 30.0;

    vector<vector<index_t>> delta_regions =
      normal_delta_regions_greedy(inds, part, normal_degrees_delta);
    vector<surface> surfaces;
    for (auto r : delta_regions) {
      surfaces.push_back(surface(&part, r));
    }
    return surfaces;
  }

  bool
  share_edge(const std::vector<gca::edge>& edges,
	     const surface& l,
	     const surface& r) {
    vector<gca::edge> l_es = edges;
    delete_if(l_es,
	      [l](const gca::edge e) { return !(l.contains(e.l) && l.contains(e.r)); });

    vector<gca::edge> r_es = edges;
    delete_if(r_es,
	      [r](const gca::edge e) { return !(r.contains(e.l) && r.contains(e.r)); });

    return intersection(l_es, r_es).size() > 0;
  }

  std::vector<surface_group>
  convex_surface_groups(const std::vector<surface*>& surfaces) {
    if (surfaces.size() == 0) { return {}; }

    vector<gca::edge> conv_edges = convex_edges(surfaces.front()->get_parent_mesh());
    cout << "# of convex edges = " << conv_edges.size() << endl;
    auto ccs =
      connected_components_by(surfaces,
			      [&conv_edges](const surface* l, const surface* r)
			      { return share_edge(conv_edges, *l, *r); });
    return ccs;
  }

  std::vector<surface>
  constant_orientation_subsurfaces(const surface& surf) {
    vector<index_t> inds = surf.index_list();
    std::vector<std::vector<index_t> > regions =
      normal_delta_regions(inds, surf.get_parent_mesh(), 3.0);
    std::vector<surface> surfs;
    for (auto r : regions) {
      surfs.push_back(surface(&surf.get_parent_mesh(), r));
    }
    return surfs;
  }

  std::vector<surface>
  inds_to_surfaces(const std::vector<std::vector<index_t>> regions,
		   const triangular_mesh& m) {
    vector<surface> sfs;
    for (auto r : regions) {
      sfs.push_back(surface(&m, r));
    }
    return sfs;
  }

  std::vector<std::vector<index_t>>
  surfaces_to_inds(const std::vector<surface>& surfs) {
    std::vector<std::vector<index_t>> inds;
    for (auto s : surfs) {
      inds.push_back(s.index_list());
    }
    return inds;
  }

  std::vector<surface>
  connected_vertical_surfaces(std::vector<index_t>& inds,
			      const triangular_mesh& m,
			      const point n) {
    delete_if(inds, [m, n](const index_t i)
	      { return !within_eps(angle_between(n, m.face_orientation(i)), 90, 2.0); });
    vector<vector<index_t>> regions = connect_regions(inds, m);
    return inds_to_surfaces(regions, m);
  }
  
  std::vector<surface>
  connected_vertical_surfaces(const triangular_mesh& m, const point n) {
    vector<index_t> inds = m.face_indexes();
    return connected_vertical_surfaces(inds, m, n);
  }


  boost::optional<surface>
  mesh_top_surface(const triangular_mesh& m, const point n) {
    auto surfs = outer_surfaces(m);

    cout << "Outer surfaces" << endl;
    //vtk_debug_highlight_inds(surfs);

    delete_if(surfs,
	      [n](const surface& s) { return !s.parallel_to(n, 1.0); });

    cout << "Base surfaces = " << endl;

    //vtk_debug_highlight_inds(surfs);

    if (surfs.size() > 0) {
      return merge_surfaces(surfs);
    } else {
      return boost::none;
    }
  }

  std::vector<gca::edge>
  shared_edges(const surface& r, const surface& l) {
    return intersection(r.edges(), l.edges());
  }


  point normal(const surface& s) {
    return s.face_orientation(s.front());
  }

  std::vector<gca::edge>
  boundary_edges(const surface& s) {
    std::vector<gca::edge> bound_edges;
    for (auto e : s.get_parent_mesh().edges()) {
      auto l_face_neighbors = s.get_parent_mesh().vertex_face_neighbors(e.l);
      auto r_face_neighbors = s.get_parent_mesh().vertex_face_neighbors(e.r);
      auto face_neighbors = intersection(l_face_neighbors, r_face_neighbors);
      bool contains_some_neighbors = false;
      bool contains_all_neighbors = true;
      for (auto facet : face_neighbors) {
	if (s.contains(facet)) {
	  contains_some_neighbors = true;
	} else {
	  contains_all_neighbors = false;
	}
      }
      if (contains_some_neighbors && !contains_all_neighbors) {
	bound_edges.push_back(e);
      }
    }
    return bound_edges;
  }

  plane surface_plane(const surface& s) {
    triangle t = s.face_triangle(s.front());
    return plane(t.normal, t.v1);
  }

  surface find_surface_by_normal(const std::vector<surface>& surfs,
				 const point n) {
    auto r = find_if(begin(surfs), end(surfs),
		     [n](const surface& s)
		     { return within_eps(angle_between(normal(s), n), 0, 1.0); });
    DBG_ASSERT(r != end(surfs));
    return *r;
  }


  std::vector<gca::edge>
  orthogonal_boundary_edges(const surface& s,
			    const point n) {
    vector<gca::edge> edges;
    for (auto e : boundary_edges(s)) {
      point ev = s.vertex(e.l) - s.vertex(e.r);
      if (within_eps(angle_between(ev, n), 90.0, 1.0)) {
	edges.push_back(e);
      }
    }
    return edges;
  }

  double max_in_dir(const surface& mesh,
		    const point dir) {
    return max_distance_along(mesh.vertex_list(), dir);
  }

  double min_in_dir(const surface& mesh,
		    const point dir) {
    return min_distance_along(mesh.vertex_list(), dir);
  }

  point max_point_in_dir(const surface& mesh,
			  const point dir) {
    return max_along(mesh.vertex_list(), dir);
  }

  point min_point_in_dir(const surface& mesh,
			  const point dir) {
    return min_along(mesh.vertex_list(), dir);
  }

  bool share_orthogonal_valley_edge(const surface& l, const surface& r) {
    vector<shared_edge> shared =
      all_shared_edges(l.index_list(), r.index_list(), l.get_parent_mesh());
    for (auto s : shared) {
      if (is_valley_edge(s, l.get_parent_mesh()) &&
	  angle_eps(s, l.get_parent_mesh(), 90.0, 0.5)) {
	return true;
      }
    }

    return false;
  }

  std::vector<plane> max_area_basis(const std::vector<surface>& surfaces) {
    DBG_ASSERT(surfaces.size() > 0);
    
    vector<surface> sorted_part_surfaces = surfaces;
    sort(begin(sorted_part_surfaces), end(sorted_part_surfaces),
	 [](const surface& l, const surface& r)
	 { return l.surface_area() > r.surface_area(); });
    
    vector<surface> basis =
      take_basis(sorted_part_surfaces,
		 [](const surface& l, const surface& r)
		 { return within_eps(angle_between(normal(l), normal(r)), 90, 2.0); },
		 2);

    vector<plane> planes;
    for (auto s : basis) {
      cout << "surface normal = " << normal(s) << endl;
      //vtk_debug_highlight_inds(s);
      planes.push_back(plane(normal(s), s.face_triangle(s.front()).v1));
    }

    point third_vec = cross(planes[0].normal(), planes[1].normal());
    const triangular_mesh& m = surfaces.front().get_parent_mesh();
    point third_pt = max_point_in_dir(m, third_vec);

    planes.push_back(plane(third_vec, third_pt));

    return planes;
  }

  bool share_non_fully_concave_edge(const surface& l, const surface& r) {
    vector<shared_edge> shared =
      all_shared_edges(l.index_list(), r.index_list(), l.get_parent_mesh());

    for (auto s : shared) {
      if (is_valley_edge(s, l.get_parent_mesh())) {
	return true;
      } else if (angle_between_normals(s, l.get_parent_mesh()) < 70.0) {
	return true;
      }
    }

    return false;
  }
  
}
