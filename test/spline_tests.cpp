#include "utils/arena_allocator.h"
#include "catch.hpp"
#include "geometry/b_spline.h"

namespace gca {

  TEST_CASE("Degree 1 spline test") {
    arena_allocator a;
    set_system_allocator(&a);

    vector<point> control_points;
    control_points.push_back(point(1, 0, 0));
    control_points.push_back(point(2, 0, 0));
      
    vector<double> knots;
    knots.push_back(0.0);
    knots.push_back(0.0);
    knots.push_back(1.0);
    knots.push_back(1.0);

    int degree = 1;

    b_spline s(degree, control_points, knots);

    SECTION("basis N0,0(0.0)") {
      REQUIRE(s.basis(0, 0, 0) == 0);
    }

    SECTION("start point") {
      point res = s.eval(0.0);
      REQUIRE(within_eps(res, point(1, 0, 0)));
    }

    SECTION("Spline point 0.25") {
      REQUIRE(within_eps(s.eval(0.25), point(1.25, 0, 0)));
    }

    SECTION("Spline point 0.5") {
      REQUIRE(within_eps(s.eval(0.5), point(1.5, 0, 0)));
    }

    SECTION("end point") {
      point res = s.eval(0.999999999);
      REQUIRE(within_eps(res, point(2, 0, 0)));
    }
    
  }

}
