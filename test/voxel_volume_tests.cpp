#include "catch.hpp"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkIntArray.h>
#include <vtkDataArray.h>
#include <vtkOBBTree.h>

#include "geometry/voxel_volume.h"
#include "geometry/voxel_volume_debug.h"
#include "geometry/vtk_debug.h"
#include "geometry/vtk_utils.h"
#include "system/parse_stl.h"

namespace gca {

  point min_point(const box b) {
    return point(b.x_min, b.y_min, b.z_min);
  }

  bool contains_point_obb(const point pt,
			  vtkSmartPointer<vtkOBBTree>& tree) {
 
    // Intersect the locator with the line
    double lineP0[3] = {pt.x, pt.y, pt.z};
    double lineP1[3] = {-10000.0, -1000.0, -1000.0};
    double lineP2[3] = {100.0, 100.0, 100.0};

    vtkSmartPointer<vtkPoints> intersectPoints = 
      vtkSmartPointer<vtkPoints>::New();
 
    int res = tree->IntersectWithLine(lineP0, lineP1, intersectPoints, NULL);
    if (res != -1) {
      return false;
    }

    int res2 = tree->IntersectWithLine(lineP0, lineP2, intersectPoints, NULL);

    if (res == -1 && res2 == -1) {
      return true;
    } else {
      //cout << "res = " << res << endl;
      //assert(res == 1);
      return false;
    }

 
    // std::cout << "NumPoints: " << intersectPoints->GetNumberOfPoints()
    // 	      << std::endl;
 
    // // Display list of intersections
    // double intersection[3];
    // for(int i = 0; i < intersectPoints->GetNumberOfPoints(); i++ )
    //   {
    // 	intersectPoints->GetPoint(i, intersection);
    // 	std::cout << "Intersection " << i << ": "
    // 		  << intersection[0] << ", "
    // 		  << intersection[1] << ", "
    // 		  << intersection[2] << std::endl;
    //   }
     
  }

  int num_line_intersections(const point pt0,
			     const point pt1,
			     vtkSmartPointer<vtkOBBTree>& tree) {
 
    // Intersect the locator with the line
    double lineP0[3] = {pt0.x, pt0.y, pt0.z};
    double lineP1[3] = {pt1.x, pt1.y, pt1.z};

    vtkSmartPointer<vtkPoints> intersectPoints = 
      vtkSmartPointer<vtkPoints>::New();
    int res2 = tree->IntersectWithLine(lineP0, lineP1, intersectPoints, NULL);

    return intersectPoints->GetNumberOfPoints();
    // std::cout << "NumPoints: " << intersectPoints->GetNumberOfPoints()
    // 	      << std::endl;
    
    // if (res == -1 && res2 == -1) {
    //   return true;
    // } else {
    //   //cout << "res = " << res << endl;
    //   //assert(res == 1);
    //   return false;
    // }

  }
  
  bool contains_point(const point p,
		      vtkSmartPointer<vtkPolyData>& pd) {
    double testInside[3] = {p.x, p.y, p.z};

    vtkSmartPointer<vtkPoints> points = 
      vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(testInside);
 
    vtkSmartPointer<vtkPolyData> pointsPolydata = 
      vtkSmartPointer<vtkPolyData>::New();
    pointsPolydata->SetPoints(points);
 
    //Points inside test
    vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = 
      vtkSmartPointer<vtkSelectEnclosedPoints>::New();
    selectEnclosedPoints->SetInputData(pointsPolydata);
    selectEnclosedPoints->SetSurfaceData(pd);
    selectEnclosedPoints->Update();
 
    return selectEnclosedPoints->IsInside(0);
  }

  voxel_volume build_from_mesh(const triangular_mesh& m) {
    box bb = m.bounding_box();
    double resolution = bb.x_len() / 80.0;

    auto pd = polydata_for_trimesh(m);
    
    DBG_ASSERT(is_closed(pd));

    vtkSmartPointer<vtkOBBTree> tree = 
      vtkSmartPointer<vtkOBBTree>::New();
    tree->SetDataSet(pd);
    tree->BuildLocator();
    
    voxel_volume vv(min_point(bb), bb.x_len(), bb.y_len(), bb.z_len(), resolution);

    for (int i = 0; i < vv.num_x_elems(); i++) {
      for (int j = 0; j < vv.num_y_elems(); j++) {
	for (int k = 0; k < vv.num_z_elems(); k++) {

	  point pt(vv.x_center(i), vv.y_center(j), vv.z_center(k));
	  if (contains_point_obb(pt, tree)) {
	    vv.set_occupied(i, j, k);
	  }

	}
      }
    }
    
    return vv;
  }

  voxel_volume accessible_from_direction(const triangular_mesh& m,
					 const point access_dir) {
    auto pd = polydata_for_trimesh(m);
    
    DBG_ASSERT(is_closed(pd));

    vtkSmartPointer<vtkOBBTree> tree = 
      vtkSmartPointer<vtkOBBTree>::New();
    tree->SetDataSet(pd);
    tree->BuildLocator();

    box bb = m.bounding_box();
    double resolution = bb.x_len() / 40.0;
    voxel_volume vv(min_point(bb), bb.x_len(), bb.y_len(), bb.z_len(), resolution);

    for (int i = 0; i < vv.num_x_elems(); i++) {
      for (int j = 0; j < vv.num_y_elems(); j++) {
	for (int k = 0; k < vv.num_z_elems(); k++) {

	  point pt(vv.x_center(i), vv.y_center(j), vv.z_center(k));
	  point outside_mesh = pt + (bb.x_len()*20)*access_dir;

	  if (num_line_intersections(pt, outside_mesh, tree) == 0) {
	    vv.set_occupied(i, j, k);
	  }

	}
      }
    }
    
    return vv;
  }

  voxel_volume difference(const voxel_volume& top,
			  const voxel_volume& bottom) {

    DBG_ASSERT(top.num_x_elems() == bottom.num_x_elems());
    DBG_ASSERT(top.num_y_elems() == bottom.num_y_elems());
    DBG_ASSERT(top.num_z_elems() == bottom.num_z_elems());
    
    voxel_volume diff(top.get_origin(),
		      top.x_length(),
		      top.y_length(),
		      top.z_length(),
		      top.get_resolution());

    
    for (int i = 0; i < top.num_x_elems(); i++) {
      for (int j = 0; j < top.num_y_elems(); j++) {
	for (int k = 0; k < top.num_z_elems(); k++) {
	  if (top.is_occupied(i, j, k) &&
	      !bottom.is_occupied(i, j, k)) {
	    diff.set_occupied(i, j, k);
	  }
	}
      }
    }

    return diff;
  }

  TEST_CASE("Initial voxel volume is empty") {
    voxel_volume vol(point(0, 0, 0), 1.0, 1.0, 1.0, 0.1);

    REQUIRE(vol.is_empty(0, 0, 0));
  }


  TEST_CASE("After setting the voxel it is occupied") {
    voxel_volume vol(point(0, 0, 0), 1.0, 1.0, 1.0, 0.1);

    vol.set_occupied(0, 0, 0);

    REQUIRE(vol.is_occupied(0, 0, 0));

  }

    vector<string> voxel_test_parts{
      "test/stl-files/onshape_parts/Part Studio 1 - Part 1(3).stl",
	"test/stl-files/onshape_parts/Part Studio 1 - Part 1(20).stl",
	"test/stl-files/onshape_parts/Part Studio 1 - Part 1(24).stl", 
	"test/stl-files/onshape_parts/Part Studio 1 - Part 1(2).stl", 
	"test/stl-files/onshape_parts/PSU Mount - PSU Mount.stl",
	"test/stl-files/onshape_parts/Part Studio 1 - Part 1(29).stl", 
	"test/stl-files/OctagonWithHolesShort.stl",
	"test/stl-files/CircleWithFilletAndSide.stl",
	"test/stl-files/onshape_parts/100-013 - Part 1.stl",
	"test/stl-files/onshape_parts/Part Studio 1 - ESC spacer.stl",
	"test/stl-files/onshape_parts/Part Studio 1 - Part 1(23).stl",
	"test/stl-files/onshape_parts/Japanese_Two_Contours_Part.stl",
	"test/stl-files/onshape_parts/Part Studio 1 - Part 1.stl",
	"test/stl-files/onshape_parts/Part Studio 1 - Falcon Prarie .177 single shot tray.stl"};
  
  TEST_CASE("Loading a model from an stl file") {
    triangular_mesh m =
      //parse_stl("test/stl-files/onshape_parts/PSU Mount - PSU Mount.stl", 0.0001);
      //parse_stl("test/stl-files/onshape_parts/Magnetic Latch Top - Part 1.stl", 0.0001);
      parse_stl("test/stl-files/onshape_parts/Part Studio 1 - Part 1(10).stl",
		0.0001);

    for (auto& n : voxel_test_parts) {
      triangular_mesh m = parse_stl(n, 0.0001);
      vtk_debug_mesh(m);

      voxel_volume bottom = accessible_from_direction(m, point(0, 1, 0));

      cout << "done with point 0 0 -1" << endl;
      vtk_debug_voxel_volume(bottom);

      voxel_volume top = accessible_from_direction(m, point(1, 0, 0));

      vtk_debug_voxel_volume(top);
      //build_from_mesh(m);

      voxel_volume diff = difference(top, bottom);

      vtk_debug_voxel_volume(diff);
    }
  }

}
