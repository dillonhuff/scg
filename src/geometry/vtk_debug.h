#ifndef GCA_VTK_DEBUG_H
#define GCA_VTK_DEBUG_H

#include <vtkActor.h>
#include <vtkPlane.h>
#include <vtkPolyData.h>
#include <vtkPlaneSource.h>
#include <vtkSmartPointer.h>

#include "geometry/depth_field.h"
#include "geometry/rigid_arrangement.h"
#include "geometry/surface.h"

namespace gca {

  class color {
  protected:
    unsigned r, g, b;

  public:
    color(unsigned p_r, unsigned p_g, unsigned p_b) :
      r(p_r), g(p_g), b(p_b) {}

    unsigned red() const { return r; }
    unsigned green() const { return g; }
    unsigned blue() const { return b; }
  };

  void highlight_cells(vtkSmartPointer<vtkPolyData> polyData,
		       const std::vector<std::pair<std::vector<index_t >, color> >& inds);
  
  color random_color(const color mix);

  bool has_cell_normals(vtkPolyData* polydata);
  bool is_closed(vtkPolyData* polydata);
  
  void debug_print_summary(vtkPolyData* polydata);
  void debug_print_edge_summary(vtkPolyData* pdata);
  void debug_print_is_closed(vtkPolyData* polydata);
  void debug_print_polydata(vtkPolyData* polydata);

  vtkSmartPointer<vtkActor> polydata_actor(vtkSmartPointer<vtkPolyData> polyData);
  void visualize_actors(const std::vector<vtkSmartPointer<vtkActor> >& pds);
  void visualize_polydatas(const std::vector<vtkSmartPointer<vtkPolyData> >& actors);

  void vtk_debug_triangles(const std::vector<triangle>& mesh);
  void vtk_debug_highlight_inds(const std::vector<index_t>& inds,
				const triangular_mesh& mesh);
  void vtk_debug_mesh(const triangular_mesh& mesh);

  void vtk_debug_meshes(const std::vector<const triangular_mesh*>& mesh);

  void vtk_debug_meshes(const std::vector<triangular_mesh>& meshes);

  void vtk_debug_highlight_inds(const surface& surf);
  void vtk_debug_highlight_inds(const std::vector<surface>& surfs);

  void vtk_debug_polygon(const oriented_polygon& p);
  void vtk_debug_polygons(const std::vector<oriented_polygon>& polys);

  void vtk_debug_mesh_boundary_edges(const triangular_mesh& m);
  void debug_arrangement(const rigid_arrangement& a);
  void vtk_debug(const triangular_mesh& m, const plane pl);

  vtkSmartPointer<vtkActor> plane_actor(vtkSmartPointer<vtkPlane> pl);

  void vtk_debug_ring(const std::vector<point>& pts);

  vtkSmartPointer<vtkActor> plane_actor(vtkSmartPointer<vtkPlane> pl);

  void vtk_debug_polygon(const labeled_polygon_3& p);
  void vtk_debug_polygons(const std::vector<labeled_polygon_3>& p);

  std::vector<vtkSmartPointer<vtkActor>>
  polygon_3_actors(const polygon_3& p);

  vtkSmartPointer<vtkActor>
  mesh_actor(const triangular_mesh& m);

  vtkSmartPointer<vtkPolyData>
  polydata_for_depth_field(const depth_field& df);
  
  void vtk_debug_depth_field(const depth_field& df);

  color random_color_non_pastel(const color mix);

  void
  visualize_surface_decomp(const std::vector<std::vector<surface> >& surf_complexes);
  
}

#endif
