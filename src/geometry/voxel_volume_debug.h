#pragma once

#include "geometry/voxel_volume.h"
#include "geometry/vtk_debug.h"

namespace gca {

  void vtk_debug_voxel_volume(const voxel_volume& df);

  vtkSmartPointer<vtkPolyData>
  polydata_for_voxel_volume(const voxel_volume& vv);
  
}
