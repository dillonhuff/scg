#ifndef GCA_PARSE_STL_H
#define GCA_PARSE_STL_H

#include "geometry/triangle.h"
#include "geometry/triangular_mesh.h"

namespace gca {

  struct stl_data {
    string name;
    vector<triangle> triangles;

    stl_data(string namep) : name(namep) {}
  };

  stl_data parse_stl(const string& stl_path);

  triangular_mesh parse_stl(const string& stl_path, const double tolerance);

  triangular_mesh parse_and_scale_stl(const std::string& part_path,
				      const double scale_factor,
				      const double tol);

  triangular_mesh parse_and_scale_box_stl(const std::string& part_path,
					  const double max_dim,
					  const double tol);
  
}

#endif
