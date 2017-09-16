#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkFeatureEdges.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkClipPolyData.h>
#include <vtkProperty.h>
#include <vtkCellArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkVertex.h>

#include <vtkAxesActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>

#include "geometry/vtk_debug.h"
#include "geometry/vtk_utils.h"

namespace gca {

  void highlight_cells(vtkSmartPointer<vtkPolyData> polyData,
		       const std::vector<index_t>& inds) {
    vtkSmartPointer<vtkUnsignedCharArray> colors = 
      vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");
 
    for(index_t i = 0; i < polyData->GetNumberOfCells(); i++) {
      unsigned char color[3];
      if (elem(i, inds)) {
	color[0] = 200;
	color[1] = 0;
	color[2] = 0;
      } else {
	color[0] = 0;
	color[1] = 0;
	color[2] = 200;
	
      }

      colors->InsertNextTupleValue(color);
    }
 
    polyData->GetCellData()->SetScalars(colors);
  }

  void highlight_cells(vtkSmartPointer<vtkPolyData> polyData,
		       const std::vector<std::pair<std::vector<index_t >, color> >& inds) {
    vtkSmartPointer<vtkUnsignedCharArray> colors = 
      vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");
 
    for(index_t i = 0; i < polyData->GetNumberOfCells(); i++) {
      unsigned char color[3];

      bool not_elem = true;
      for (auto& ind_color_pair : inds) {
	if (elem(i, ind_color_pair.first)) {
	  auto& cl = ind_color_pair.second;

	  color[0] = cl.red();
	  color[1] = cl.green();
	  color[2] = cl.blue();
	  not_elem = false;
	  break;
	}
      }
	
      if (not_elem) {
	color[0] = 255;
	color[1] = 255;
	color[2] = 255;
      }
      
      colors->InsertNextTupleValue(color);
    }
 
    polyData->GetCellData()->SetScalars(colors);
  }
  
  vtkSmartPointer<vtkActor> polydata_actor(vtkSmartPointer<vtkPolyData> polyData)
  {
    vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);
 
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    return actor;
  }

  vtkSmartPointer<vtkActor> plane_actor(vtkSmartPointer<vtkPlane> pl)
  {
    double n[3];
    pl->GetNormal(n);
    double pt[3];
    pl->GetOrigin(pt);
    
    vtkSmartPointer<vtkPlaneSource> planeSource =
      vtkSmartPointer<vtkPlaneSource>::New();
    planeSource->SetCenter(pt);
    planeSource->SetNormal(n);
    planeSource->Update();
 
    vtkPolyData* plane = planeSource->GetOutput();
 
    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(plane);
 
    vtkSmartPointer<vtkActor> actor =
      vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);    

    return actor;
  }

  void
  visualize_polydatas(const std::vector<vtkSmartPointer<vtkPolyData> >& pds) {
    std::vector<vtkSmartPointer<vtkActor> > actors;

    for (auto& pd : pds) {
      actors.push_back(polydata_actor(pd));
    }

    visualize_actors(actors);
  }
  
  void visualize_actors(const std::vector<vtkSmartPointer<vtkActor> >& actors)
  {
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    // Create axes
    vtkSmartPointer<vtkAxesActor> axes =
      vtkSmartPointer<vtkAxesActor>::New();

    renderer->AddActor(axes);
    for (auto actor : actors) {
      renderer->AddActor(actor);
    }
    renderer->SetBackground(0.0, 0.0, 0.0); //1.0, 1.0, 1.0); //0.0, 0.0, 0.0); //.1, .2, .3);

    vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderWindow->Render();
    renderWindowInteractor->Start();
  }

  vtkSmartPointer<vtkActor>
  highlighted_surface_actor(const std::vector<index_t>& inds,
			    const triangular_mesh& mesh) {
    auto poly_data = polydata_for_trimesh(mesh);
    highlight_cells(poly_data, inds);
    return polydata_actor(poly_data);
  }

  void vtk_debug_highlight_inds(const std::vector<index_t>& inds,
				const triangular_mesh& mesh) {
    auto surface_act = highlighted_surface_actor(inds, mesh);
    visualize_actors({surface_act});
  }

  void vtk_debug_triangles(const std::vector<triangle>& tris) {
    auto pd = polydata_for_triangles(tris);
    vtkSmartPointer<vtkActor> act = polydata_actor(pd);
    visualize_actors({act});
  }
  
  void vtk_debug_mesh(const triangular_mesh& mesh) {
    vtk_debug_highlight_inds(mesh.face_indexes(), mesh);
  }

  void vtk_debug_polygon(const oriented_polygon& p) {
    auto pd = polydata_for_polygon(p);
    debug_print_edge_summary(pd);
    auto act = polydata_actor(pd);
    visualize_actors({act});
  }

  void vtk_debug_meshes(const std::vector<const triangular_mesh*>& meshes) {
    std::vector<vtkSmartPointer<vtkActor>> actors;
    for (auto m : meshes) {
      auto p = polydata_for_trimesh(*m);
      actors.push_back(polydata_actor(p));
    }
    visualize_actors(actors);
  }

  void vtk_debug_meshes(const std::vector<triangular_mesh>& meshes) {
    std::vector<vtkSmartPointer<vtkActor>> actors;
    for (auto m : meshes) {
      auto p = polydata_for_trimesh(m);
      actors.push_back(polydata_actor(p));
    }
    visualize_actors(actors);
  }
  
  void vtk_debug_highlight_inds(const surface& surf) {
    vtk_debug_highlight_inds(surf.index_list(), surf.get_parent_mesh());
  }

  void vtk_debug_highlight_inds(const std::vector<surface>& surfs) {
    vtk_debug_highlight_inds(merge_surfaces(surfs));
  }

  bool is_closed(vtkPolyData* polydata)
  {
    vtkSmartPointer<vtkFeatureEdges> featureEdges = 
      vtkSmartPointer<vtkFeatureEdges>::New();
    featureEdges->FeatureEdgesOff();
    featureEdges->BoundaryEdgesOn();
    featureEdges->NonManifoldEdgesOn();
    featureEdges->SetInputData(polydata);
    featureEdges->Update();
 
    int number_of_open_edges =
      featureEdges->GetOutput()->GetNumberOfCells();

    return number_of_open_edges == 0;
  }

  bool has_cell_normals(vtkPolyData* polydata)
  {
    // std::cout << "Looking for cell normals..." << std::endl;
 
    // Count points
    //vtkIdType numCells = polydata->GetNumberOfCells();
    //std::cout << "There are " << numCells << " cells." << std::endl;
 
    // Count triangles
    //vtkIdType numPolys = polydata->GetNumberOfPolys();
    //std::cout << "There are " << numPolys << " polys." << std::endl;
 
    ////////////////////////////////////////////////////////////////
    // Double normals in an array
    vtkDoubleArray* normalDataDouble =
      vtkDoubleArray::SafeDownCast(polydata->GetCellData()->GetArray("Normals"));
 
    if(normalDataDouble)
      {
	//int nc = normalDataDouble->GetNumberOfTuples();
	// std::cout << "There are " << nc
	// 	  << " components in normalDataDouble" << std::endl;
	return true;
      }
 
    ////////////////////////////////////////////////////////////////
    // Double normals in an array
    vtkFloatArray* normalDataFloat =
      vtkFloatArray::SafeDownCast(polydata->GetCellData()->GetArray("Normals"));
 
    if(normalDataFloat)
      {
	//int nc = normalDataFloat->GetNumberOfTuples();
	// std::cout << "There are " << nc
	// 	  << " components in normalDataFloat" << std::endl;
	return true;
      }
 
    ////////////////////////////////////////////////////////////////
    // Point normals
    vtkDoubleArray* normalsDouble =
      vtkDoubleArray::SafeDownCast(polydata->GetCellData()->GetNormals());
 
    if(normalsDouble)
      {
	// std::cout << "There are " << normalsDouble->GetNumberOfComponents()
	// 	  << " components in normalsDouble" << std::endl;
	return true;
      }
 
    ////////////////////////////////////////////////////////////////
    // Point normals
    vtkFloatArray* normalsFloat =
      vtkFloatArray::SafeDownCast(polydata->GetCellData()->GetNormals());
 
    if(normalsFloat)
      {
	// std::cout << "There are " << normalsFloat->GetNumberOfComponents()
	// 	  << " components in normalsFloat" << std::endl;
	return true;
      }
 
    /////////////////////////////////////////////////////////////////////
    // Generic type point normals
    vtkDataArray* normalsGeneric = polydata->GetCellData()->GetNormals(); //works
    if(normalsGeneric)
      {
	std::cout << "There are " << normalsGeneric->GetNumberOfTuples()
		  << " normals in normalsGeneric" << std::endl;
 
	double testDouble[3];
	normalsGeneric->GetTuple(0, testDouble);
 
	std::cout << "Double: " << testDouble[0] << " "
		  << testDouble[1] << " " << testDouble[2] << std::endl;
 
	return true;
      }
 
    return false;
  }  
  
  void debug_print_summary(vtkPolyData* polydata) {
    cout << "# of points in polydata = " << polydata->GetNumberOfPoints() << endl;
    cout << "# of polys in polydata = " << polydata->GetNumberOfPolys() << endl;
    bool has_normals = has_cell_normals(polydata);
    cout << "Has normals ? " << (has_normals == 1 ? "True" : "False") << endl;
  }

  void debug_print_edge_summary(vtkPolyData* pdata) {
    vtkSmartPointer<vtkFeatureEdges> boundEdges = 
      vtkSmartPointer<vtkFeatureEdges>::New();
    boundEdges->FeatureEdgesOff();
    boundEdges->BoundaryEdgesOn();
    boundEdges->NonManifoldEdgesOff();
    boundEdges->SetInputData(pdata);
    boundEdges->Update();
 
    int number_of_boundary_edges =
      boundEdges->GetOutput()->GetNumberOfCells();

    vtkSmartPointer<vtkFeatureEdges> nonManifoldEdges = 
      vtkSmartPointer<vtkFeatureEdges>::New();
    nonManifoldEdges->FeatureEdgesOff();
    nonManifoldEdges->BoundaryEdgesOff();
    nonManifoldEdges->ManifoldEdgesOff();
    nonManifoldEdges->NonManifoldEdgesOn();
    nonManifoldEdges->SetInputData(pdata);
    nonManifoldEdges->Update();

    int number_of_non_manifold_edges =
      nonManifoldEdges->GetOutput()->GetNumberOfCells();

    cout << "# of boundary edges = " << number_of_boundary_edges << endl;
    cout << "# of non manifold edges = " << number_of_non_manifold_edges << endl;

    vtkIndent* indent = vtkIndent::New();
    if (number_of_non_manifold_edges > 0) {
      for (int i = 0; i < nonManifoldEdges->GetOutput()->GetNumberOfCells(); i++) {
	vtkCell* c = nonManifoldEdges->GetOutput()->GetCell(i);
	assert(c->GetCellType() == VTK_LINE);
	c->PrintSelf(cout, *indent);
      }
    }
  }

  void debug_print_is_closed(vtkPolyData* polydata) {
    if(is_closed(polydata)) {
      std::cout << "Surface is closed" << std::endl;
    } else {
      std::cout << "Surface is not closed" << std::endl;
    }
  }
  
  void debug_print_polydata(vtkPolyData* polydata) {
    debug_print_summary(polydata);
    debug_print_is_closed(polydata);

    vtkIndent* indent = vtkIndent::New();
    for (vtkIdType i = 0; i < polydata->GetNumberOfPolys(); i++) {
      vtkCell* c = polydata->GetCell(i);
      c->PrintSelf(cout, *indent);
      cout << "Cell type = " << c->GetCellType() << endl;
      cout << "# of points = " << c->GetNumberOfPoints() << endl;
    }
  }

  void vtk_debug_polygons(const std::vector<oriented_polygon>& polys) {
    vector<vtkSmartPointer<vtkActor>> actors;
    for (auto p : polys) {
      auto pd = polydata_for_polygon(p);
      actors.push_back(polydata_actor(pd));
    }
    visualize_actors(actors);
  }

  void vtk_debug_mesh_boundary_edges(const triangular_mesh& m) {
    auto pd = polydata_for_trimesh(m);

    vtkSmartPointer<vtkFeatureEdges> featureEdges = 
      vtkSmartPointer<vtkFeatureEdges>::New();
    featureEdges->FeatureEdgesOff();
    featureEdges->BoundaryEdgesOn();
    featureEdges->NonManifoldEdgesOff();
    featureEdges->SetInputData(pd);
    featureEdges->Update();

    auto edge_actor = polydata_actor(featureEdges->GetOutput());
    auto pd_actor = polydata_actor(pd);

    visualize_actors({pd_actor, edge_actor});
  }

  void debug_arrangement(const rigid_arrangement& a) {
    cout << "Debugging" << endl;
    vector<vtkSmartPointer<vtkActor>> actors;
    for (auto name : a.mesh_names()) {
      if (a.metadata(name).display_during_debugging) {
	auto pd = polydata_for_trimesh(a.mesh(name));
	actors.push_back(polydata_actor(pd));
      }
    }
    visualize_actors(actors);
  }

  void vtk_debug(const triangular_mesh& m,
		 const plane pl) {
    auto mpd = polydata_for_trimesh(m);
    auto mpl = vtk_plane(pl);

    visualize_actors({polydata_actor(mpd), plane_actor(mpl)});
  }

  vtkSmartPointer<vtkActor>
  mesh_actor(const triangular_mesh& m) {
    auto mpd = polydata_for_trimesh(m);
    return polydata_actor(mpd);
  }


  void vtk_debug_ring(const std::vector<point>& pts) {
    auto pd = polydata_for_ring(pts);
    vector<vtkSmartPointer<vtkActor>> actors{polydata_actor(pd)};
    visualize_actors(actors);
  }

  std::vector<vtkSmartPointer<vtkActor>>
  polygon_3_actors(const polygon_3& p) {
    vector<vtkSmartPointer<vtkPolyData>> ring_pds;

    auto pd = polydata_for_ring(p.vertices());
    ring_pds.push_back(pd);
    cout << "??? # lines in poly = " << pd->GetNumberOfPolys() << endl;
    
    for (auto ir : p.holes()) {
      auto pd = polydata_for_ring(ir);
      cout << "??? # lines in hole? = " << pd->GetNumberOfPolys() << endl;
      ring_pds.push_back(polydata_for_ring(ir));
    }
    
    vector<vtkSmartPointer<vtkActor>> ring_acts;
    for (auto r : ring_pds) {
      ring_acts.push_back(polydata_actor(r));
    }

    return ring_acts;
  }
  
  void vtk_debug_polygon(const labeled_polygon_3& p) {
    // vector<vtkSmartPointer<vtkPolyData>> ring_pds;

    // auto pd = polydata_for_ring(p.vertices());
    // ring_pds.push_back(pd);
    // cout << "??? # lines in poly = " << pd->GetNumberOfPolys() << endl;
    
    // for (auto ir : p.holes()) {
    //   auto pd = polydata_for_ring(ir);
    //   cout << "??? # lines in hole? = " << pd->GetNumberOfPolys() << endl;
    //   ring_pds.push_back(polydata_for_ring(ir));
    // }
    
    vector<vtkSmartPointer<vtkActor>> ring_acts = polygon_3_actors(p);
    // for (auto r : ring_pds) {
    //   ring_acts.push_back(polydata_actor(r));
    // }

    visualize_actors(ring_acts);
  }

  void vtk_debug_polygons(const std::vector<labeled_polygon_3>& polys) {
    vector<vtkSmartPointer<vtkPolyData>> ring_pds;

    for (auto p : polys) {
      auto pd = polydata_for_ring(p.vertices());
      ring_pds.push_back(pd);
    
      for (auto ir : p.holes()) {
	auto pd = polydata_for_ring(ir);
	ring_pds.push_back(polydata_for_ring(ir));
      }
    }
    
    vector<vtkSmartPointer<vtkActor>> ring_acts;
    for (auto r : ring_pds) {
      ring_acts.push_back(polydata_actor(r));
    }

    visualize_actors(ring_acts);
  }

  color random_color(const color mix) {
    unsigned red = rand() % 256;
    unsigned green = rand() % 256;
    unsigned blue = rand() % 256;

    // mix the color
    red = (red + mix.red()) / 2;
    green = (green + mix.green()) / 2;
    blue = (blue + mix.blue()) / 2;

    color color(red, green, blue);
    return color;    
  }


  void vtk_debug_depth_field(const depth_field& df) {
    auto pd_actor = polydata_actor(polydata_for_depth_field(df));

    visualize_actors({pd_actor});
  }

  vtkSmartPointer<vtkPolyData>
  polydata_for_depth_field(const depth_field& df) {
    vtkSmartPointer<vtkPoints> points =
      vtkSmartPointer<vtkPoints>::New();
 
    vtkSmartPointer<vtkCellArray> vertices =
      vtkSmartPointer<vtkCellArray>::New();
    
    for (int i = 0; i < df.num_x_elems; i++) {
      for (int j = 0; j < df.num_y_elems; j++) {

	auto id = points->InsertNextPoint(df.x_center(i),
					  df.y_center(j),
					  df.column_height(i, j));
	vtkSmartPointer<vtkVertex> vertex = 
	  vtkSmartPointer<vtkVertex>::New();

	vertex->GetPointIds()->SetId(0, id);
	vertices->InsertNextCell(vertex);
      }
    }
 
    vtkSmartPointer<vtkPolyData> polydata =
      vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetVerts(vertices);

    return polydata;
  }


  color random_color_non_pastel(const color mix) {
    unsigned red = rand() % 256;
    unsigned green = rand() % 256;
    unsigned blue = rand() % 256;

    return color(red, green, blue);

  }

  void
  visualize_surface_decomp(const std::vector<std::vector<surface> >& surf_complexes)
  {

    cout << "# of surfaces = " << surf_complexes.size() << endl;
    if (surf_complexes.size() == 0) { return; }

    const auto& part = surf_complexes.front().front().get_parent_mesh();

    auto part_polydata = polydata_for_trimesh(part);

    color white(255, 255, 255);

    vector<pair<vector<index_t>,  color> > colors;
    for (auto& sc : surf_complexes) {

      vector<index_t> inds;
      for (auto& s : sc) {
	concat(inds, s.index_list());
      }

      color tp_color = random_color_non_pastel(white);
      colors.push_back(std::make_pair(inds, tp_color));

    }

    highlight_cells(part_polydata, colors);
    visualize_actors({polydata_actor(part_polydata)});

  }

}
