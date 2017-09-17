#include "catch.hpp"
#include "geometry/plane.h"

namespace gca {

  TEST_CASE("Plane line intersection") {
    plane p(point(0, 0, 1), point(0, 0, 0));

    SECTION("Line that does not intersect") {
      line l(point(10, 0, 6), point(11, 0, 7));
      REQUIRE(!intersection(p, l));
    }

    SECTION("Line that does intersect") {
      line l(point(0, 0, 10), point(0, 0, -10));
      boost::optional<point> r = intersection(p, l);
      REQUIRE(r);
      REQUIRE(within_eps(*r, point(0, 0, 0), 0.001));
    }
  }
}
