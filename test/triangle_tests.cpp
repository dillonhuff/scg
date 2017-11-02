#include "catch.hpp"
#include "geometry/triangle.h"
#include "utils/algorithm.h"
#include "utils/arena_allocator.h"
#include "system/parse_stl.h"

namespace gca {

  TEST_CASE("Segment intersection with triangles") {
    auto triangles = parse_stl("/Users/dillon/CppWorkspace/gca/test/stl-files/SlicedCone.stl").triangles;

    SECTION("Line at 10.0") {
      point s(0.332573, 0.612317, 10.0);
      point e(0.357573, 0.112317, 10.0);
      line l(s, e);
      REQUIRE(!intersects_triangles(l, triangles));
    }

    SECTION("Line at 0.1") {
      point s(0.332573, 0.612317, 0.1);
      point e(0.357573, 0.112317, 0.1);
      line l(s, e);
      REQUIRE(intersects_triangles(l, triangles));
    }
  }
  
}
