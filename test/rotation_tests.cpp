#include "catch.hpp"
#include "geometry/rotation.h"

namespace gca {

  TEST_CASE("Simple 90 degree rotation") {
    point from(0, 1, 0);
    point to(0, 0, 1);

    const rotation r = rotate_from_to(from, to);

    REQUIRE(within_eps(determinant(r), 1.0, 0.001));
  }

  TEST_CASE("Simple 2 component 90 degree rotation") {
    point from(1, 1, 0);
    point to(0, 0, 1);

    const rotation r = rotate_from_to(from, to);

    REQUIRE(within_eps(determinant(r), 1.0, 0.001));
  }
  
  TEST_CASE("All positive rotation, not 90 degrees") {
    point from(7.59, 3.4, 9.32);
    point to(0, 0, 1);

    const rotation r = rotate_from_to(from, to);

    REQUIRE(within_eps(determinant(r), 1.0, 0.001));
  }

  TEST_CASE("All positive rotation, not 90 degrees, close to (0, 0, 1)") {
    point from(0.59, 0.4, 9.32);
    point to(0, 0, 1);

    const rotation r = rotate_from_to(from, to);

    REQUIRE(within_eps(determinant(r), 1.0, 0.001));
  }

  TEST_CASE("All positive rotation") {
    point from(0.90, 0.32, 0);
    point to(0, 0, 1);

    const rotation r = rotate_from_to(from, to);

    REQUIRE(within_eps(determinant(r), 1.0, 0.001));
  }
  
  TEST_CASE("Random rotation") {
    point from(0.90, -0.32, 0);
    point to(0, 0, 1);

    const rotation r = rotate_from_to(from, to);

    REQUIRE(within_eps(determinant(r), 1.0, 0.001));
  }
  
  TEST_CASE("Irregular 90 degree rotation") {
    point from(0.953399, -0.301714, -0);
    point to(0, 0, 1);

    const rotation r = rotate_from_to(from, to);

    REQUIRE(within_eps(determinant(r), 1.0, 0.001));
  }
  
}
