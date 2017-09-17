#include "catch.hpp"
#include "geometry/arc.h"
#include "utils/arena_allocator.h"

namespace gca {

  TEST_CASE("Arc locations") {
    arena_allocator a;
    set_system_allocator(&a);

    SECTION("Start location") {
      arc a = arc(point(1, 0, 0), point(1, 1, 0), point(0, 0.5, 0), CLOCKWISE);
      REQUIRE(within_eps(a.value(0.0), point(1, 0, 0)));
    }

    SECTION("Vertical arc") {
      arc a = arc(point(2, 1, 1), point(2, 2, 1), point(0, 0.5, 0), CLOCKWISE);
      SECTION("Start") {
	REQUIRE(within_eps(a.value(0.0), point(2, 1, 1)));
      }

      SECTION("End") {
	REQUIRE(within_eps(a.value(1.0), point(2, 2, 1)));
      }

      SECTION("Middle") {
	REQUIRE(within_eps(a.value(0.5), point(1.5, 1.5, 1)));
      }
      
    }

    SECTION("Horizontal arc") {
      arc a = arc(point(0, 0, 0), point(1, 0, 0), point(0.5, 0, 0), COUNTERCLOCKWISE);
      SECTION("Start") {
	REQUIRE(within_eps(a.value(0.0), point(0, 0, 0)));
      }

      SECTION("End") {
	REQUIRE(within_eps(a.value(1.0), point(1, 0, 0)));
      }

      SECTION("Middle") {
	REQUIRE(within_eps(a.value(0.5), point(0.5, -0.5, 0)));
      }
    }

    SECTION("Random example originally from gcode_to_cuts") {
      arc a(point(1, 1, 1),
	    point(3, 4.5, 1),
	    point(1, 1.75, 0),
	    COUNTERCLOCKWISE);

      SECTION("Start") {
	point s = a.value(0.0);
	REQUIRE(within_eps(s, point(1, 1, 1)));
      }

      SECTION("End") {
	REQUIRE(within_eps(a.value(1.0), point(3, 4.5, 1)));
      }

      SECTION("Middle") {
	point m = a.value(0.5);
	REQUIRE(within_eps(m, point(3.75, 1.75, 1)));
      }
    }

  }
}
