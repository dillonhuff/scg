#include "catch.hpp"
#include "geometry/line.h"
#include "geometry/polyline.h"
#include "utils/arena_allocator.h"

namespace gca {

  TEST_CASE("Offset line creation") {
    arena_allocator a;
    set_system_allocator(&a);
    
    SECTION("One segment") {
      polyline p({point(0, 0, 0), point(1, 0, 0)});
      auto off = offset(p, OFFSET_LEFT, 3.0);
      polyline correct({point(0, 3, 0), point(1, 3, 0)});
      REQUIRE(pointwise_within_eps(off, correct, 0.000001));
    }

    SECTION("Unconnected convex and concave segment") {
      point p1(0, 0, 0);
      point p2(0, 1, 0);
      point p3(2, 1, 0);
      point p4(3, 2, 0);
      double inc = 0.5;
      double deg = -90;
      polyline p({p1, p2, p3, p4});
      auto off = offset(p, OFFSET_RIGHT, inc);
      point c1(0.5, 0, 0);
      point c2(0.5, 0.5, 0);
      line l1(p2, p3);
      line l2(p3, p4);
      point i = inc*((p3 - p2).normalize().rotate_z(deg));      
      line k1(l1.start + i, l1.end + i);
      i = inc*((p4 - p3).normalize().rotate_z(deg));
      line k2(l2.start + i, l2.end + i);
      point c3 = trim_or_extend_unsafe(k1, k2);
      point c4 = p4 + i;
      polyline correct({c1, c2, c3, c4});
      REQUIRE(pointwise_within_eps(off, correct, 0.00001));
    }

    SECTION("Fully connected triangles") {
      point p1(1, 0, -1);
      point p2(4, -23, -1);
      point p3(10, 3, -1);
      polyline p({p1, p2, p3, p1});

      double inc = 0.2;
      double deg = 90;

      auto off = offset(p, OFFSET_LEFT, inc);

      line l1(p1, p2);
      point i1 = inc*((p2 - p1).normalize().rotate_z(deg));
      line l2(p2, p3);
      point i2 = inc*((p3 - p2).normalize().rotate_z(deg));
      line l3(p3, p1);
      point i3 = inc*((p1 - p3).normalize().rotate_z(deg));
      
      line k1 = l1.shift(i1);
      line k2 = l2.shift(i2);
      line k3 = l3.shift(i3);

      point c1 = trim_or_extend_unsafe(k3, k1);
      point c2 = trim_or_extend_unsafe(k1, k2);
      point c3 = trim_or_extend_unsafe(k2, k3);

      polyline correct({c1, c2, c3, c1});
      REQUIRE(pointwise_within_eps(off, correct, 0.00001));
    }

    SECTION("Offset that produced nan") {
      point p1(0, 0, 0);
      point p2(1, 0, 0);
      point p3(1, -1, 0);

      polyline p({p1, p2, p3});
      double inc = 0.1;
      double deg = 90;

      auto off = offset(p, OFFSET_LEFT, inc);

      point c1(0, 0.1, 0);
      point c2(1.1, 0.1, 0);
      point c3(1.1, -1, 0);
      polyline correct({c1, c2, c3});

      cout << "Offset square" << endl;
      for (auto pt : off) {
      	cout << pt << endl;
      }

      REQUIRE(pointwise_within_eps(off, correct, 0.00001));
    }

    SECTION("Adjacent parallel segments") {
      point p1(0, 0, 3.4);
      point p2(0, 1, 2.1);
      point p3(0, 2, 0);

      polyline p({p1, p2, p3});
      double inc = 1.0;
      double deg = 90;

      auto off = offset(p, OFFSET_LEFT, inc);

      point c1(-1, 0, 3.4);
      point c2(-1, 1, 2.1);
      point c3(-1, 2, 0);
      polyline correct({c1, c2, c3});

      cout << "Offset square" << endl;
      for (auto pt : off) {
      	cout << pt << endl;
      }

      REQUIRE(pointwise_within_eps(off, correct, 0.00001));
    }
  }

  TEST_CASE("Exterior direction") {
    polyline p({point(0, 0, 0),
	  point(1, 0, 0),
	  point(1, 1, 0),
	  point(0, 1, 0),
	  point(0, 0, 0)});

    SECTION("exterior direction is OFFSET_RIGHT") {
      REQUIRE(exterior_direction(p) == OFFSET_RIGHT);
    }

    SECTION("exterior direction is OFFSET_LEFT") {
      reverse(begin(p), end(p));
      REQUIRE(exterior_direction(p) == OFFSET_LEFT);
    }

  }
}
