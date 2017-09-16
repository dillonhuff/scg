#ifndef GCA_EXTRUSION_H
#define GCA_EXTRUSION_H

#include "geometry/polygon_3.h"
#include "geometry/vtk_utils.h"

namespace gca {

  typedef std::vector<index_t> index_poly;

  struct extrusion {
    std::vector<point> pts;
    std::vector<index_poly> poly_layers;
    std::vector<double> layer_depths;
    point extrude_dir;
  };

  std::vector<index_poly>
  unordered_segments_to_index_polygons(std::vector<gca::edge>& lines);

  std::vector<index_poly>
  unordered_segments_to_index_polylines(std::vector<gca::edge>& lines);
  
  triangular_mesh
  extrude_layers(const std::vector<point>& pts,
		 const std::vector<index_poly>& poly_layers,
		 const std::vector<double>& layer_depths,
		 const point extrude_dir);

  triangular_mesh
  extrude(const extrusion& ext);

  polyline
  to_polyline(const index_poly& poly,
	      const std::vector<point>& pts);

  oriented_polygon
  oriented_polygon_for_index_polyline(const std::vector<point>& pts,
				      const index_poly& p,
				      const point n);

  triangular_mesh extrude(const polygon_3& p, const point v);

  triangular_mesh
  extrude_surface_negative(const std::vector<index_t>& surf,
			   const triangular_mesh& part,
			   const point n,
			   const double length);

}

#endif
