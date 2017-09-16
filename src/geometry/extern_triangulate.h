#pragma once

#include <vector>

#include "geometry/polygon_3.h"
#include "geometry/triangle.h"

namespace gca {

  std::vector<triangle> triangulate_flat_3d(const polygon_3& p);

}
