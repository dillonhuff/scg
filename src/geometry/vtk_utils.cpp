#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolygon.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>

#include "geometry/vtk_debug.h"
#include "geometry/vtk_utils.h"

namespace gca {

  line vtkCell_to_line(vtkCell* c) {
    DBG_ASSERT(c->GetCellType() == VTK_LINE);
    DBG_ASSERT(c->GetNumberOfPoints() == 2);

    vtkLine* ln = dynamic_cast<vtkLine*>(c);
    double p0[3];
    double p1[3];
    ln->GetPoints()->GetPoint(0, p0);
    ln->GetPoints()->GetPoint(1, p1);

    point v0(p0[0], p0[1], p0[2]);
    point v1(p1[0], p1[1], p1[2]);

    return line(v0, v1);
  }

  triangle vtkCell_to_triangle(vtkCell* c) {
    DBG_ASSERT(c->GetCellType() == VTK_TRIANGLE);
    DBG_ASSERT(c->GetNumberOfPoints() == 3);

    vtkTriangle* tri = dynamic_cast<vtkTriangle*>(c);
    double p0[3];
    double p1[3];
    double p2[3];
    tri->GetPoints()->GetPoint(0, p0);
    tri->GetPoints()->GetPoint(1, p1);
    tri->GetPoints()->GetPoint(2, p2);

    point v0(p0[0], p0[1], p0[2]);
    point v1(p1[0], p1[1], p1[2]);
    point v2(p2[0], p2[1], p2[2]);

    point q1 = v1 - v0;
    point q2 = v2 - v0;
    point norm = cross(q2, q1).normalize();
    return triangle(norm, v0, v1, v2);
  }

  std::vector<triangle>
  polydata_to_triangle_list(vtkPolyData* in_polydata) {
    vtkSmartPointer<vtkTriangleFilter> triangleFilter =
      vtkSmartPointer<vtkTriangleFilter>::New();
    triangleFilter->SetInputData(in_polydata);
    triangleFilter->Update();

    vtkPolyData* polydata = triangleFilter->GetOutput();

    vector<triangle> tris;
    for (vtkIdType i = 0; i < polydata->GetNumberOfPolys(); i++) {
      vtkCell* c = polydata->GetCell(i);
      tris.push_back(vtkCell_to_triangle(c));
    }
    return tris;
  }

  vtkSmartPointer<vtkPolyData>
  polydata_for_triangles(const std::vector<triangle>& tris) {
    vtkSmartPointer<vtkPoints> points =
      vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> triangles =
      vtkSmartPointer<vtkCellArray>::New();

    for (auto t : tris) {
      vtkIdType i1 = points->InsertNextPoint(t.v1.x, t.v1.y, t.v1.z);
      vtkIdType i2 = points->InsertNextPoint(t.v2.x, t.v2.y, t.v2.z);
      vtkIdType i3 = points->InsertNextPoint(t.v3.x, t.v3.y, t.v3.z);

      vtkSmartPointer<vtkTriangle> triangle =
	vtkSmartPointer<vtkTriangle>::New();

      triangle->GetPointIds()->SetId(0, i1);
      triangle->GetPointIds()->SetId(1, i2);
      triangle->GetPointIds()->SetId(2, i3);

      triangles->InsertNextCell(triangle);    
    }

    // Create a polydata object
    vtkSmartPointer<vtkPolyData> polyData =
      vtkSmartPointer<vtkPolyData>::New();
 
    // Add the geometry and topology to the polydata
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);

    return polyData;
  }

  void color_polydata(vtkSmartPointer<vtkPolyData> data,
		      const unsigned char red,
		      const unsigned char green,
		      const unsigned char blue) {
    unsigned char color[3];
    color[0] = red;
    color[1] = green;
    color[2] = blue;

    // Create a vtkUnsignedCharArray container and store the colors in it
    vtkSmartPointer<vtkUnsignedCharArray> colors =
      vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    for (vtkIdType i = 0; i < data->GetNumberOfCells(); i++) {
      colors->InsertNextTupleValue(color);
    }
 

    data->GetCellData()->SetScalars(colors);

  }
  
  vtkSmartPointer<vtkPolyData>
  polydata_for_ring(const std::vector<point>& mpts) {
    // Create the polydata where we will store all the geometric data
    vtkSmartPointer<vtkPolyData> linesPolyData =
      vtkSmartPointer<vtkPolyData>::New();
 
 
    // Create a vtkPoints container and store the points in it
    vtkSmartPointer<vtkPoints> pts =
      vtkSmartPointer<vtkPoints>::New();
    for (auto p : mpts) {
      pts->InsertNextPoint(p.x, p.y, p.z);
    }
 
    // Add the points to the polydata container
    linesPolyData->SetPoints(pts);
 
    vtkSmartPointer<vtkCellArray> lines =
      vtkSmartPointer<vtkCellArray>::New();

    for (vtkIdType i = 0; i < linesPolyData->GetNumberOfPoints(); i++) { 
      vtkSmartPointer<vtkLine> line0 =
	vtkSmartPointer<vtkLine>::New();
      line0->GetPointIds()->SetId(0, i);
      line0->GetPointIds()->SetId(1, (i + 1) % linesPolyData->GetNumberOfPoints());

      lines->InsertNextCell(line0);
    }
 
    // Add the lines to the polydata container
    linesPolyData->SetLines(lines);
 
 
    unsigned char red[3] = { 255, 0, 0 };
 
    // Create a vtkUnsignedCharArray container and store the colors in it
    vtkSmartPointer<vtkUnsignedCharArray> colors =
      vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    for (vtkIdType i = 0; i < lines->GetNumberOfCells(); i++) {
      colors->InsertNextTupleValue(red);
    }
 

    linesPolyData->GetCellData()->SetScalars(colors);

    return linesPolyData;
  }
  
  vtkSmartPointer<vtkPolyData>
  polydata_for_trimesh(const triangular_mesh& mesh) {
    vtkSmartPointer<vtkPoints> points =
      vtkSmartPointer<vtkPoints>::New();

    // TODO: Does i++ matter here?
    for (auto i : mesh.vertex_indexes()) {
      point p = mesh.vertex(i);
      points->InsertNextPoint(p.x, p.y, p.z);
    }

    vtkSmartPointer<vtkCellArray> triangles =
      vtkSmartPointer<vtkCellArray>::New();

    for (auto i : mesh.face_indexes()) {
      vtkSmartPointer<vtkTriangle> triangle =
	vtkSmartPointer<vtkTriangle>::New();

      auto t = mesh.triangle_vertices(i);
      triangle->GetPointIds()->SetId(0, t.v[0]);
      triangle->GetPointIds()->SetId(1, t.v[1]);
      triangle->GetPointIds()->SetId(2, t.v[2]);

      triangles->InsertNextCell(triangle);    
    }

    // Create a polydata object
    vtkSmartPointer<vtkPolyData> polyData =
      vtkSmartPointer<vtkPolyData>::New();
 
    // Add the geometry and topology to the polydata
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);

    return polyData;
  }

  triangular_mesh
  trimesh_for_polydata(vtkPolyData* in_polydata) {
    auto tris = polydata_to_triangle_list(in_polydata);

    vtk_debug_triangles(tris);

    triangular_mesh m = make_mesh(tris, 0.00001);

    DBG_ASSERT(m.is_connected());

    return m;
  }

  vtkSmartPointer<vtkPolyData>
  polydata_for_polygon(const oriented_polygon& p) {
    vtkSmartPointer<vtkPoints> points =
      vtkSmartPointer<vtkPoints>::New();

    for (auto vert : p.vertices()) {
      points->InsertNextPoint(vert.x, vert.y, vert.z);
    }

    vtkSmartPointer<vtkPolygon> pl = vtkSmartPointer<vtkPolygon>::New();
    pl->GetPointIds()->SetNumberOfIds(p.vertices().size());
    for (int i = 0; i < p.vertices().size(); i++) {
      if (!(pl->GetPointIds())) {
	cout << "No point ids in pl!" << endl;
	DBG_ASSERT(false);
      }
      pl->GetPointIds()->SetId(i, i);
    }

    vtkSmartPointer<vtkCellArray> polys =
      vtkSmartPointer<vtkCellArray>::New();
    polys->InsertNextCell(pl);
    
    // Create a polydata object
    vtkSmartPointer<vtkPolyData> polyData =
      vtkSmartPointer<vtkPolyData>::New();
 
    // Add the geometry and topology to the polydata
    polyData->SetPoints(points);
    polyData->SetPolys(polys);

    return polyData;
  }

  std::vector<triangle>
  vtk_triangulate_poly(const oriented_polygon& p) {
    auto polyData =
      polydata_for_polygon(p);

    return polydata_to_triangle_list(polyData);

  }

  vtkSmartPointer<vtkPlane> vtk_plane(const plane pl) {
    vtkSmartPointer<vtkPlane> clipPlane = 
      vtkSmartPointer<vtkPlane>::New();
    point n = pl.normal();
    point pt = pl.pt();
    clipPlane->SetNormal(n.x, n.y, n.z);
    clipPlane->SetOrigin(pt.x, pt.y, pt.z);
    return clipPlane;
  }

  void append_polyline(unsigned npts,
		       vtkSmartPointer<vtkPolyData> linesPolyData,
		       vtkSmartPointer<vtkPoints> pts,
		       vtkSmartPointer<vtkCellArray> lines,
		       const polyline& pl) {
    for (auto p : pl) 
      {
	pts->InsertNextPoint(p.x, p.y, p.z);
      }

    if (pl.num_points() > 0) 
      {
	
    	for (vtkIdType i = 0; i < pl.num_points() - 1; i++) { 
    	  vtkSmartPointer<vtkLine> line0 =
    	    vtkSmartPointer<vtkLine>::New();
    	  line0->GetPointIds()->SetId(0, npts + i);
    	  line0->GetPointIds()->SetId(1, npts + i + 1);

    	  lines->InsertNextCell(line0);
    	}
      }
    
 
  }

  vtkSmartPointer<vtkPolyData>
  polydata_for_polylines(const std::vector<polyline>& polylines) {
    // Create the polydata where we will store all the geometric data
    vtkSmartPointer<vtkPolyData> linesPolyData =
      vtkSmartPointer<vtkPolyData>::New();
 
    // Create a vtkPoints container and store the points in it
    vtkSmartPointer<vtkPoints> pts =
      vtkSmartPointer<vtkPoints>::New();

    // Add the points to the polydata container
    linesPolyData->SetPoints(pts);
    vtkSmartPointer<vtkCellArray> lines =
      vtkSmartPointer<vtkCellArray>::New();

    // Add the lines to the polydata container
    linesPolyData->SetLines(lines);

    for (auto& pl : polylines) {
      auto num_points_so_far = linesPolyData->GetNumberOfPoints();
      append_polyline(num_points_so_far, linesPolyData, pts, lines, pl);
    }

    cout << "# of lines = " << linesPolyData->GetNumberOfLines() << endl;

    return linesPolyData;
  }
  
}
