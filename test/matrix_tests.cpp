#include "catch.hpp"
#include "geometry/matrix.h"

namespace gca {

  TEST_CASE("Create rotation matrix") {
    point a(0, 0, 1);
    point b(1, 0, 0);
    point c(0, 1, 0);

    point ap(0, -1, 0);
    point bp(1, 0, 0);
    point cp(0, 0, 1);

    const ublas::matrix<double> r =
      plane_basis_rotation(a, b, c,
			   ap, bp, cp);

    REQUIRE(within_eps(determinant(r), 1.0));

    REQUIRE(within_eps(times_3(r, a), ap));
    REQUIRE(within_eps(times_3(r, b), bp));
    REQUIRE(within_eps(times_3(r, c), cp));
  }
}
