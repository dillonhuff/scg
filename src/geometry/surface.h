#ifndef GCA_SURFACE_H
#define GCA_SURFACE_H

#include <boost/optional.hpp>

#include "geometry/plane.h"
#include "geometry/triangular_mesh.h"

namespace gca {

  class surface {
  protected:
    const triangular_mesh* parent_mesh;
    vector<index_t> tri_indexes;

  public:
    inline bool contains(index_t ind) const {
      return std::binary_search(begin(tri_indexes), end(tri_indexes), ind);
    }

    inline bool contains_vertex(index_t ind) const {
      auto neighbor_faces = get_parent_mesh().vertex_face_neighbors(ind);
      for (auto i : neighbor_faces) {
	if (contains(i)) { return true; }
      }
      return false;
    }

    point vertex(const index_t i) const { return get_parent_mesh().vertex(i); }

    bool is_boundary_edge(const gca::edge e) const {
      auto face_neighbs = get_parent_mesh().edge_face_neighbors(e);
      for (auto n : face_neighbs) {
	if (!contains(n)) { return true; }
      }
      return false;
    }

    std::vector<index_t>
    edge_face_neighbors(const gca::edge e) const {
      auto m_neighbs = get_parent_mesh().edge_face_neighbors(e);
      vector<index_t> neighbs;
      for (auto i : m_neighbs) {
	if (contains(i)) {
	  neighbs.push_back(i);
	}
      }
      return neighbs;
    }
    
    std::vector<gca::edge> edges() const {
      std::vector<gca::edge> edgs;
      for (auto e : get_parent_mesh().edges()) {
	if (this->edge_face_neighbors(e).size() > 0) {
	  edgs.push_back(e);
	}
      }
      return edgs;
    }

    std::vector<index_t> index_list() const {
      return tri_indexes;
    }

    std::vector<point> vertex_list() const {
      vector<point> verts;
      for (auto i : index_list()) {
	triangle t = get_parent_mesh().face_triangle(i);
	verts.push_back(t.v1);
	verts.push_back(t.v2);
	verts.push_back(t.v3);
	//	verts.push_back(vertex(i));
      }

      return verts;
    }

    bool contained_by_sorted(const std::vector<index_t>& inds) const {
      for (auto i : tri_indexes) {
	if (!binary_search(begin(inds), end(inds), i)) {
	  return false;
	}
      }
      return true;
    }

    bool contained_by(const surface& other) const {
      return contained_by_sorted(other.tri_indexes);
    }
    
    inline index_t front() const {
      return tri_indexes.front();
    }

    double surface_area() const {
      double total = 0.0;
      for (auto i : tri_indexes) {
	total += (parent_mesh->face_triangle(i)).area();
      }
      return total;
    }

    inline point face_orientation(index_t ind) const
    { return parent_mesh->face_orientation(ind); }

    inline triangle face_triangle(index_t ind) const
    { return get_parent_mesh().face_triangle(ind); }
    
    surface(const triangular_mesh* p_parent_mesh,
	    const vector<index_t>& p_tri_indexes) :
      parent_mesh(p_parent_mesh), tri_indexes(p_tri_indexes) {
      assert(tri_indexes.size() > 0);
      std::sort(begin(tri_indexes), end(tri_indexes));
    }

    inline const triangular_mesh& get_parent_mesh() const
    { return *parent_mesh; }

    bool orthogonal_to(const point n, double tol) const {
      return gca::all_orthogonal_to(index_list(), get_parent_mesh(), n, tol);
    }

    bool parallel_to(const point n, double tol) const {
      return gca::all_parallel_to(index_list(), get_parent_mesh(), n, tol);
    }

    bool antiparallel_to(const point n, double tol) const {
      return gca::all_antiparallel_to(index_list(), get_parent_mesh(), n, tol);
    }
    
  };

  bool surfaces_share_edge(const unsigned i,
			   const unsigned j,
			   const std::vector<surface>& surfaces);

  bool surfaces_share_edge(const unsigned i,
			   const unsigned j,
			   const std::vector<surface*>& surfaces);
  
  bool surfaces_share_edge(const surface& l,
			   const surface& r);

  
  std::vector<index_t> surface_vertexes(const surface& s);

  bool orthogonal_flat_surfaces(const surface* l, const surface* r);
  bool parallel_flat_surfaces(const surface* l, const surface* r);
  std::vector<surface> outer_surfaces(const triangular_mesh& part);

  typedef std::vector<std::vector<index_t>> surface_list;

  surface merge_surfaces(const std::vector<surface>& surfaces);

  void remove_contained_surfaces(const std::vector<surface>& stable_surfaces,
				 std::vector<surface>& surfaces_to_cut);

  boost::optional<oriented_polygon>
  part_outline(std::vector<surface>* surfaces_to_cut);

  boost::optional<surface>
  part_outline_surface(std::vector<surface>* surfaces_to_cut,
		       const point n);

  typedef std::vector<unsigned> surface_group;

  boost::optional<surface>
  part_outline_surface(const triangular_mesh& m,
		       const point n);

  std::vector<surface> surfaces_to_cut(const triangular_mesh& part);

  bool
  share_edge(const std::vector<gca::edge>& edges,
	     const surface& l,
	     const surface& r);

  std::vector<surface_group>
  convex_surface_groups(const std::vector<surface*>& surfaces);

  std::vector<surface>
  constant_orientation_subsurfaces(const surface& surface);

  std::vector<surface>
  connected_vertical_surfaces(const triangular_mesh& m, const point n);

  std::vector<surface>
  connected_vertical_surfaces(std::vector<index_t>& inds,
			      const triangular_mesh& m,
			      const point n);


  boost::optional<surface>
  mesh_top_surface(const triangular_mesh& m, const point n);

  std::vector<surface> surfaces_to_cut(const std::vector<index_t>& indexes,
				       const triangular_mesh& part);

  std::vector<surface>
  merge_surface_groups(const std::vector<surface>& surfs,
		       const std::vector<std::vector<unsigned> >& groups);

  std::vector<gca::edge>
  shared_edges(const surface& r, const surface& l);

  std::vector<gca::edge>
  boundary_edges(const surface& r);
  
  point normal(const surface& s);

  plane surface_plane(const surface& s);

  surface find_surface_by_normal(const std::vector<surface>& surfs,
				 const point n);

  std::vector<std::vector<index_t>>
  surfaces_to_inds(const std::vector<surface>& surfs);

  std::vector<surface>
  inds_to_surfaces(const std::vector<std::vector<index_t>> regions,
		   const triangular_mesh& m);
  
  std::vector<gca::edge>
  orthogonal_boundary_edges(const surface& s,
			    const point n);

  double min_in_dir(const surface& mesh, const point dir);
  double max_in_dir(const surface& mesh, const point dir);

  point min_point_in_dir(const surface& mesh, const point dir);
  point max_point_in_dir(const surface& mesh, const point dir);

  bool share_orthogonal_valley_edge(const surface& l, const surface& r);

  std::vector<plane> max_area_basis(const std::vector<surface>& surfaces);

  bool share_non_fully_concave_edge(const surface& l, const surface& r);

}

#endif
