#pragma once

#include "geometry/polygon.h"
#include "geometry/polygon_3.h"
#include "geometry/triangular_mesh.h"

namespace gca {

  struct shared_edge {
    index_t triangle_1;
    index_t triangle_2;
    edge e;
  };

  bool angle_eps(const shared_edge e,
		 const triangular_mesh& m,
		 const double angle,
		 const double tol);

  bool is_valley_edge(const shared_edge e,
		      const triangular_mesh& m);

  std::vector<shared_edge> all_shared_edges(const std::vector<index_t>& l_faces,
					    const std::vector<index_t>& r_faces,
					    const triangular_mesh& part);

  vector<oriented_polygon> mesh_bounds(const vector<index_t>& faces,
				       const triangular_mesh& mesh);

  vector<polygon_3> surface_boundary_polygons(const vector<index_t>& faces,
					      const triangular_mesh& mesh);
  
  oriented_polygon max_area_outline(const std::vector<index_t>& inds,
				    const triangular_mesh& m);

  polygon_3 box_bound(const triangular_mesh& mesh);

  double angle_between_normals(const shared_edge e,
			       const triangular_mesh& m);

  polygon_3 surface_boundary_polygon(const vector<index_t>& faces,
				     const triangular_mesh& mesh);

  std::vector<index_t>
  vertex_inds_on_surface(const std::vector<index_t>& s,
			 const triangular_mesh& m);
}
