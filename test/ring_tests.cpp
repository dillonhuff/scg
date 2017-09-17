#include "catch.hpp"
#include "geometry/ring.h"

namespace gca {

  TEST_CASE("Ring clockwise direction") {

    vector<point> safe_ring{point(-3, -3, 0),
	point(3, -3, 0),
	point(3, 3, 0),
	point(-3, 3, 0)};

    point clockwise_dir(0, 0, -1);

    point clock_normal = clockwise_normal(safe_ring);

    cout << "Computed clockwise normal = " << clock_normal << endl;
    REQUIRE(within_eps(clock_normal, clockwise_dir));
  }

  TEST_CASE("Ring counterclockwise direction") {

    vector<point> safe_ring{point(-3, -3, 0),
	point(3, -3, 0),
	point(3, 3, 0),
	point(-3, 3, 0)};

    point clockwise_dir(0, 0, 1);

    point counterclock_normal = counterclockwise_normal(safe_ring);

    cout << "Computed counterclockwise normal = " << counterclock_normal << endl;
    REQUIRE(within_eps(counterclock_normal, clockwise_dir));
  }
  
}
