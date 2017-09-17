#include "catch.hpp"
#include "geometry/line.h"

namespace gca {

  TEST_CASE("Same segment several times") {
    point bl(0, 0, 0);
    point ul(0, 1, 0);
    point br(1, 0, 0);
    point ur(1, 1, 0);

    vector<line> lines{line(ur, bl), line(ur, bl)};
    REQUIRE(count_in(line(ur, bl), lines) == 2);
  }

  TEST_CASE("Line intersection") {

    SECTION("Parallel horizontal lines") {
      point a(0, 0, 0);
      point b(1, 0, 0);

      point c(1, 1, 0);
      point d(2, 1, 0);

      line l1(a, b);
      line l2(c, d);

      auto res = line_intersection_2d(l1, l2);
      REQUIRE(!res.just);
    }

    SECTION("Horizontal and vertical lines") {
      point a(0, 0, 0);
      point b(1, 0, 0);

      point c(2, 3, 0);
      point d(2, 2, 0);

      line l1(a, b);
      line l2(c, d);

      point correct(2, 0, 0);

      SECTION("Horizontal, vertical") {
	auto res = line_intersection_2d(l1, l2);
	REQUIRE(res.just);
	REQUIRE(within_eps(correct, res.t));
      }

      SECTION("Vertical, horizontal") {
	auto res = line_intersection_2d(l2, l1);
	REQUIRE(res.just);
	REQUIRE(within_eps(correct, res.t));
      }

    }
  }
}
