#include "catch.hpp"
#include "geometry/point.h"

namespace gca {

  TEST_CASE("Rotate XY 90 degrees counterclockwise") {
    point p(0, 1, 0);
    point r = p.rotate_z(90);
    point c(-1, 0, 0);
    REQUIRE(within_eps(r, c));
  }

  TEST_CASE("Rotate XY 90 degrees clockwise") {
    point p(1, 0, 0);
    point r = p.rotate_z(-90);
    point c(0, -1, 0);
    REQUIRE(within_eps(r, c));
  }

  TEST_CASE("Extend back horizontal vector") {
    point s(1, 0, 0);
    point e(2, 0, 0);
    point sp = extend_back(s, e, 1);
    point c = point(0, 0, 0);
    REQUIRE(within_eps(sp, c));
  }

  TEST_CASE("Extend back (1, 1)") {
    point s(0, 0, 0);
    point e(1, 1, 0);
    point sp = extend_back(s, e, 1);
    double v = sqrt(1.0/2.0);
    point c = point(-v, -v, 0);
    REQUIRE(within_eps(sp,c ));
  }

  TEST_CASE("Extend back (-1, 1)") {
    point s(0, 0, 0);
    point e(-1, 1, 0);
    point sp = extend_back(s, e, 1);
    double v = sqrt(1.0/2.0);
    point c = point(v, -v, 0);
    REQUIRE(within_eps(sp,c ));
  }

  TEST_CASE("Extend back (-1, -1)") {
    point s(0, 0, 0);
    point e(-1, -1, 0);
    point sp = extend_back(s, e, 1);
    double v = sqrt(1.0/2.0);
    point c = point(v, v, 0);
    REQUIRE(within_eps(sp,c ));
  }

  TEST_CASE("Normalize (-5, 0, 0)") {
    point s(-5, 0, 0);
    point r = s.normalize();
    point c = point(-1, 0, 0);
    REQUIRE(within_eps(c, r));
  }

  TEST_CASE("Angle between is 180") {
    point p(1, 1, 0);
    point q(-1, -1, 0);
    double a = angle_between(p, q);
    REQUIRE(within_eps(a, 180, 0.00001));
  }

  TEST_CASE("Angle between is 90") {
    point p(1, 0, 0);
    point q(0, 0, 1);
    double a = angle_between(p, q);
    REQUIRE(within_eps(a, 90, 0.00001));
  }
  
  TEST_CASE("Angle between is 180 2") {
    point p(-1, -1.75, 0);
    point q(1, 1.75, 0);
    double a = angle_between(p, q);
    cout << "angle between = " << a << endl;
    REQUIRE(within_eps(a, 180, 0.00001));
  }

  TEST_CASE("Angle between is 180 3") {
    point p(0, -1, 0);
    point q(0, 1, 0);
    double a = angle_between(p, q);
    cout << "angle between = " << a << endl;
    REQUIRE(within_eps(a, 180, 0.00001));
  }
  
  TEST_CASE("Cross product") {
    point a(1.5, 2.7, -3.4);
    point b(5.3, -2.9, 6.0);
    point c(6.34, -27.02, -18.66);
    REQUIRE(within_eps(cross(a, b), c));
  }

  TEST_CASE("Signed distance along") {
    point n(0, 0, -1);
    point q(0, 0, 0.5);
    point s(0, 0, -0.5);
    REQUIRE(signed_distance_along(q, s) < signed_distance_along(s, n));
  }

}
