#include <vtkMassProperties.h>
#include <vtkSTLWriter.h>
#include <vtkImplicitDataSet.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkDecimatePro.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkStripper.h>
#include <vtkCutter.h>
#include <vtkFeatureEdges.h>
#include <vtkProperty.h>
#include <vtkFeatureEdges.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkClipPolyData.h>
#include <vtkProperty.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataNormals.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>

#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include "geometry/mesh_operations.h"
#include "geometry/surface.h"
#include "geometry/vtk_debug.h"
#include "geometry/vtk_utils.h"

namespace gca {

  template <class HDS>
  class build_mesh : public CGAL::Modifier_base<HDS> {
  public:
    const triangular_mesh& m;
    build_mesh(const triangular_mesh& p_m) : m(p_m) {}
  
    void operator()( HDS& hds) {
      // Postcondition: hds is a valid polyhedral surface.
      CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);

      unsigned num_verts = m.vertex_indexes().size();
      unsigned num_faces = m.face_indexes().size();

      B.begin_surface( num_verts, num_faces, 6*num_faces);

      typedef typename HDS::Vertex   Vertex;
      typedef typename Vertex::Point Point;

      for (auto i : m.vertex_indexes()) {
	point p = m.vertex(i);
	B.add_vertex(Point(p.x, p.y, p.z));
      }

      for (auto i : m.face_indexes()) {
	triangle_t t = m.triangle_vertices(i);
	B.begin_facet();
	B.add_vertex_to_facet(t.v[0]);
	B.add_vertex_to_facet(t.v[1]);
	B.add_vertex_to_facet(t.v[2]);
	B.end_facet();
      }

      B.end_surface();

    }
  };

  std::vector<triangle> polyhedron_to_triangle_list(const Polyhedron& poly) {
    vector<triangle> tris;
    for (auto it = poly.facets_begin(); it != poly.facets_end(); it++) {
      auto fc = *it;

      auto he = fc.facet_begin();
      auto end = fc.facet_begin();
      vector<point> pts;
      CGAL_For_all(he, end) {
	auto h = *he;
	auto v = h.vertex();
	auto r = (*v).point();

	point res_pt(CGAL::to_double(r.x()),
		     CGAL::to_double(r.y()),
		     CGAL::to_double(r.z()));

	pts.push_back(res_pt);
      }

      DBG_ASSERT(pts.size() == 3);

      point v0 = pts[0];
      point v1 = pts[1];
      point v2 = pts[2];

      point q1 = v1 - v0;
      point q2 = v2 - v0;
      point norm = cross(q1, q2).normalize();
      tris.push_back(triangle(norm, v0, v1, v2));
    }

    return tris;
  }

  std::vector<triangle> nef_to_triangle_list(const Nef_polyhedron& p) {
    //    DBG_DBG_ASSERT(p.is_simple());

    Polyhedron poly;
    p.convert_to_polyhedron(poly);

    return polyhedron_to_triangle_list(poly);
  }

  boost::optional<triangular_mesh>
  nef_polyhedron_to_trimesh(const Nef_polyhedron& p) {

    if (!p.is_simple()) { return boost::none; }

    auto tris = nef_to_triangle_list(p);

    if (tris.size() == 0) { return boost::none; }

    return make_mesh(tris, 1e-20);
  }

  std::vector<triangular_mesh>
  nef_polyhedron_to_trimeshes(const Nef_polyhedron& p) {

    auto tris = nef_to_triangle_list(p);

    if (tris.size() == 0) { return {}; }

    return make_meshes(tris, 1e-20);
  }
  
  Nef_polyhedron trimesh_to_nef_polyhedron(const triangular_mesh& m) {
    cout << "Converting to nef" << endl;
    //vtk_debug_mesh(m);

    build_mesh<HalfedgeDS> mesh_builder(m);

    Polyhedron P;
    P.delegate(mesh_builder);

    Nef_polyhedron N(P);

    DBG_ASSERT(N.is_simple());

    cout << "Done converting to nef" << endl;

    return N.regularization();
  }

  boost::optional<std::vector<point>>
  merge_center(const std::vector<point>& l, const std::vector<point>& r) {
    if (components_within_eps(l.back(), r.front(), 0.01)) {
      std::vector<point> rest(begin(r) + 1, end(r));
      std::vector<point> lc = l;
      concat(lc, rest);
      return lc;
    }
    return boost::none;
  }
  
  boost::optional<std::vector<point>>
  merge_adjacent_chains(const std::vector<point>& l, const std::vector<point>& r) {
    auto res = merge_center(l, r);
    if (res) { return *res; }
    res = merge_center(r, l);
    if (res) { return *res; }
    std::vector<point> lc = l;
    reverse(begin(lc), end(lc));
    res = merge_center(lc, r);
    if (res) { return *res; }
    std::vector<point> rc = r;
    reverse(begin(rc), end(rc));
    res = merge_center(l, rc);
    if (res) { return *res; }
    return boost::none;
  }

  bool
  try_to_merge_chains(std::vector<std::vector<point>>& plines) {
    for (unsigned i = 0; i < plines.size(); i++) {
      std::vector<point>* pi = &(plines[i]);
      for (unsigned j = 0; j < plines.size(); j++) {
	std::vector<point>* pj = &(plines[j]);
	if (i != j) {
	  boost::optional<std::vector<point>> res = merge_adjacent_chains(*pi, *pj);
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
  
  // TODO: Templatize and merge with unordered_segments_to_index_polylines code?
  std::vector<std::vector<point>>
  connect_chains(std::vector<std::vector<point>>& chains) {
    bool merged_one = true;
    while (merged_one) {
      merged_one = try_to_merge_chains(chains);
    }
    return chains;
  }

  std::vector<oriented_polygon>
  mesh_cross_section(const triangular_mesh& m,
		     const plane p) {
    auto clipPlane = vtk_plane(p);

    vtkSmartPointer<vtkPolyData> input_mesh_data =
      polydata_for_trimesh(m);

    // Clip the source with the plane
    vtkSmartPointer<vtkClipPolyData> clipper = 
      vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetInputData(input_mesh_data);
    clipper->SetClipFunction(clipPlane);
    clipper->Update();
    
    vtkPolyData* m_data = clipper->GetOutput();

    // Now extract clipped edges
    vtkSmartPointer<vtkFeatureEdges> boundaryEdges =
      vtkSmartPointer<vtkFeatureEdges>::New();
    boundaryEdges->SetInputData(m_data);
    boundaryEdges->BoundaryEdgesOn();
    boundaryEdges->FeatureEdgesOff();
    boundaryEdges->NonManifoldEdgesOff();
    boundaryEdges->ManifoldEdgesOff();
    boundaryEdges->Update();

    vtkPolyData& edges = *(boundaryEdges->GetOutput());

    cout << "# of edges = " << edges.GetNumberOfCells() << endl;

    vector<vector<point>> ln;
    for (vtkIdType i = 0; i < edges.GetNumberOfCells(); i++) {
      vtkCell* c = edges.GetCell(i);
      line l = vtkCell_to_line(c);
      ln.push_back({l.start, l.end});
    }

    cout << "Number of initial edges = " << ln.size() << endl;

    ln = connect_chains(ln);

    cout << "Number of connected chains = " << ln.size() << endl;

    // for (auto c : ln) {
    //   cout << "CHAIN" << endl;
    //   for (auto p : c) {
    // 	cout << "---" << p << endl;
    //   }
    // }

    vector<oriented_polygon> polys;
    for (auto pt : ln) {
      auto r = pt;

      DBG_ASSERT(components_within_eps(r.front(), r.back(), 0.01));

      r.pop_back();
      
      polys.push_back(oriented_polygon(point(0, 0, 1), r));
    }

    return polys;
  }

  triangular_mesh
  clip_mesh_exact(const triangular_mesh& m,
		  const plane pl) {
    plane box_pl = pl.flip();

    point p1 = pl.pt();

    DBG_ASSERT(false);
  }
  
  // TOOD: Simplify this atrocity of a function
  triangular_mesh
  clip_mesh(const triangular_mesh& m,
	    const plane pl) {

    auto clipPlane = vtk_plane(pl);

    vtkSmartPointer<vtkPolyData> input_mesh_data =
      polydata_for_trimesh(m);

    // auto act = plane_actor(clipPlane);
    // auto pdact = polydata_actor(input_mesh_data);

    // visualize_actors({pdact, act});
    
    
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator =
      vtkSmartPointer<vtkPolyDataNormals>::New();

    normalGenerator->SetInputData(input_mesh_data);
    normalGenerator->SetSplitting(0);
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOn();
    normalGenerator->Update();

    input_mesh_data = normalGenerator->GetOutput();
    
    // Clip the source with the plane
    vtkSmartPointer<vtkClipPolyData> clipper = 
      vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetInputData(input_mesh_data);
    clipper->SetClipFunction(clipPlane);
    clipper->Update();
    
    vtkPolyData* m_data = clipper->GetOutput();

    cout << "PRE CAPPING EDGES" << endl;
    debug_print_edge_summary(m_data);
    auto m_data_act = polydata_actor(m_data);
    visualize_actors({m_data_act});

    vector<triangle> main_tris = polydata_to_triangle_list(m_data);

    auto main_mesh = make_mesh(main_tris, 0.001);
    vtk_debug_mesh(main_mesh);

    // Now extract clipped edges
    vtkSmartPointer<vtkFeatureEdges> boundaryEdges =
      vtkSmartPointer<vtkFeatureEdges>::New();
    boundaryEdges->SetInputData(m_data);
    boundaryEdges->BoundaryEdgesOn();
    boundaryEdges->FeatureEdgesOff();
    boundaryEdges->NonManifoldEdgesOff();
    boundaryEdges->ManifoldEdgesOff();
 
    vtkSmartPointer<vtkStripper> boundaryStrips =
      vtkSmartPointer<vtkStripper>::New();
    boundaryStrips->SetInputConnection(boundaryEdges->GetOutputPort());
    boundaryStrips->Update();
 
    // Change the polylines into polygons
    vtkSmartPointer<vtkPolyData> boundaryPoly =
      vtkSmartPointer<vtkPolyData>::New();
    boundaryPoly->SetPoints(boundaryStrips->GetOutput()->GetPoints());
    boundaryPoly->SetPolys(boundaryStrips->GetOutput()->GetLines());

    // Triangulate the polygons
    vtkSmartPointer<vtkTriangleFilter> triangleFilter =
      vtkSmartPointer<vtkTriangleFilter>::New();
    triangleFilter->SetInputData(boundaryPoly);
    triangleFilter->Update();

    vtkPolyData* patch_data = triangleFilter->GetOutput();

    vector<triangle> patch_tris = polydata_to_triangle_list(patch_data);
    concat(main_tris, patch_tris);

    triangular_mesh intermediate_mesh = make_mesh_no_winding_check(main_tris, 0.01);

    auto pdata = polydata_for_trimesh(intermediate_mesh);

    vtkSmartPointer<vtkPolyDataNormals> normalGenerator2 =
      vtkSmartPointer<vtkPolyDataNormals>::New();

    normalGenerator2->SetInputData(pdata);
    normalGenerator2->SetSplitting(0);
    normalGenerator2->ComputePointNormalsOn();
    normalGenerator2->ComputeCellNormalsOn();
    normalGenerator2->Update();

    pdata = normalGenerator2->GetOutput();

    debug_print_edge_summary(pdata);

    DBG_ASSERT(is_closed(pdata));

    double target =
      1.0 - (static_cast<double>(pdata->GetNumberOfPolys()) - 12.0) / 12.0;
    vtkSmartPointer<vtkDecimatePro> decimate =
      vtkSmartPointer<vtkDecimatePro>::New();
    decimate->SetInputData(pdata);
    decimate->SetTargetReduction(target);
    decimate->Update();

    pdata = decimate->GetOutput();

    vtkSmartPointer<vtkPolyDataNormals> normalGenerator3 =
      vtkSmartPointer<vtkPolyDataNormals>::New();

    normalGenerator3->SetInputData(pdata);
    normalGenerator3->SetSplitting(0);
    normalGenerator3->ComputePointNormalsOn();
    normalGenerator3->ComputeCellNormalsOn();
    normalGenerator3->Update();

    pdata = normalGenerator3->GetOutput();

    DBG_ASSERT(is_closed(pdata));
    DBG_ASSERT(has_cell_normals(pdata));

    triangular_mesh final_mesh = trimesh_for_polydata(pdata);

    DBG_ASSERT(final_mesh.is_connected());

    return final_mesh;
  }

  boost::optional<triangular_mesh>
  boolean_difference(const triangular_mesh& a, const triangular_mesh& b) {
    auto a_nef = trimesh_to_nef_polyhedron(a);
    auto b_nef = trimesh_to_nef_polyhedron(b);
    auto res = a_nef - b_nef;

    return nef_polyhedron_to_trimesh(res);
  }

  boost::optional<triangular_mesh>
  boolean_intersection(const triangular_mesh& a, const triangular_mesh& b) {
    auto a_nef = trimesh_to_nef_polyhedron(a);
    auto b_nef = trimesh_to_nef_polyhedron(b);
    auto res = a_nef.intersection(b_nef);

    return nef_polyhedron_to_trimesh(res);
  }

  std::vector<triangular_mesh>
  boolean_difference(const std::vector<triangular_mesh>& as,
		     const std::vector<triangular_mesh>& bs) {
    vector<triangular_mesh> res;

    for (auto a : as) {
      concat(res, boolean_difference(a, bs));
    }

    return res;
  }
  
  std::vector<triangular_mesh>
  boolean_difference(const triangular_mesh& a,
		     const std::vector<triangular_mesh>& bs) {
    cout << "Converting a" << endl;
    Nef_polyhedron res = trimesh_to_nef_polyhedron(a);
    cout << "Done converting a" << endl;

    box a_box = a.bounding_box();
    
    for (auto b : bs) {
      box b_box = b.bounding_box();

      if (overlap(a_box, b_box)) {
	Nef_polyhedron b_nef = trimesh_to_nef_polyhedron(b);

	res = res - b_nef;
      }
    }

    auto result_meshes = nef_polyhedron_to_trimeshes(res);

    return result_meshes;
  }

  // TODO: Delete? Not sure this is ever used
  void write_mesh_as_stl(const triangular_mesh& m,
			 const std::string& file_name) {
    auto mesh_data = polydata_for_trimesh(m);
    
    vtkSmartPointer<vtkSTLWriter> stlWriter =
      vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName(file_name.c_str());
    stlWriter->SetInputData(mesh_data);
    stlWriter->Write();
  }

  double volume(const triangular_mesh& m) {
    auto input1 = polydata_for_trimesh(m);
    
    vtkSmartPointer<vtkMassProperties> mass =
      vtkMassProperties::New();
    mass->SetInputData(input1);
    mass->Update();

    double m_volume = mass->GetVolume();

    return m_volume;
  }

  triangular_mesh nef_to_single_trimesh(const Nef_polyhedron& nef) {
    auto meshes = nef_polyhedron_to_trimeshes(nef);

    if (meshes.size() != 1) {
      cout << "# of meshes = " << meshes.size() << endl;
      for (auto m : meshes) {
	vtk_debug_mesh(m);
      }

      cout << "All meshes" << endl;
      vtk_debug_meshes(meshes);

      DBG_ASSERT(meshes.size() == 1);
    }

    return meshes.front();
  }

  triangular_mesh nef_to_single_merged_trimesh(const Nef_polyhedron& p) {
    vector<triangle> tris = nef_to_triangle_list(p);
    return make_merged_mesh(tris, 1e-20);
  }

  struct exact_volume::volume_impl {
    Nef_polyhedron nef;

    volume_impl(const Nef_polyhedron& p_nef)
      : nef(p_nef) {}

  };

  exact_volume::exact_volume(const triangular_mesh& mesh) {
    impl = new volume_impl(trimesh_to_nef_polyhedron(mesh));
  }

  exact_volume::exact_volume(volume_impl* p_impl) {
    impl = p_impl;
  }
  
  exact_volume exact_volume::subtract(const exact_volume& other) const {
    Nef_polyhedron n = impl->nef - other.impl->nef;
    volume_impl* impl = new volume_impl(n);
    return exact_volume(impl);
  }

  exact_volume exact_volume::intersection(const exact_volume& other) const {
    Nef_polyhedron n = (impl->nef).intersection(other.impl->nef);
    volume_impl* impl = new volume_impl(n);
    return exact_volume(impl);
  }

  exact_volume exact_volume::regularization() const {
    Nef_polyhedron n = (impl->nef).regularization();
    volume_impl* impl = new volume_impl(n);
    return exact_volume(impl);
  }
  
  exact_volume::~exact_volume() {
    delete impl;
  }

  std::vector<triangular_mesh> exact_volume::to_trimeshes() const {
    return nef_polyhedron_to_trimeshes(impl->nef);
  }

  triangular_mesh exact_volume::to_single_trimesh() const {
    return nef_to_single_trimesh(impl->nef);
  }

  triangular_mesh exact_volume::to_single_merged_trimesh() const {
    return nef_to_single_merged_trimesh(impl->nef);
  }


}
