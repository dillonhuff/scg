#include "catch.hpp"
#include "geometry/surface.h"
#include "utils/arena_allocator.h"
#include "system/parse_stl.h"

namespace gca {

  TEST_CASE("Arm joint has 12 vertical surfaces in (0, 1, 0)") {
    arena_allocator a;
    set_system_allocator(&a);

    triangular_mesh m = parse_stl("/Users/dillon/CppWorkspace/gca/test/stl-files/Arm_Joint_Top.stl", 0.001);

    vector<surface> verts = connected_vertical_surfaces(m, point(0, 1, 0));
    REQUIRE(verts.size() == 12);
  }
}
