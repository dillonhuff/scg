#include "catch.hpp"
#include "geometry/triangular_mesh.h"
#include "utils/arena_allocator.h"
#include "system/parse_stl.h"

namespace gca {

  TEST_CASE("Get triangles in to and out of mesh") {
    arena_allocator a;
    set_system_allocator(&a);

    SECTION("Box mesh in and out") {
      auto triangles = parse_stl("/Users/dillon/CppWorkspace/gca/test/stl-files/SlicedCone.stl").triangles;
      triangular_mesh m = make_mesh(triangles, 0.001);
      auto res_triangles = m.triangle_list();
      REQUIRE(res_triangles.size() == triangles.size());
    }
  }

  TEST_CASE("Winding order consistency test") {
    arena_allocator a;
    set_system_allocator(&a);

    auto triangles = parse_stl("/Users/dillon/CppWorkspace/gca/test/stl-files/SlicedCone.stl").triangles;
    triangular_mesh m = make_mesh(triangles, 0.001);
    REQUIRE(m.winding_order_is_consistent());
  }
}
