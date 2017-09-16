#ifndef GCA_RING_H
#define GCA_RING_H

#include "geometry/point.h"
#include "utils/check.h"

namespace gca {

  template<typename Ring>
  point ring_normal(const Ring& r) {
    DBG_ASSERT(r.size() >= 3);

    double normal_x = 0.0;
    double normal_y = 0.0;
    double normal_z = 0.0;

    int i, j;
    for (i = 0, j = 1; i < r.size(); i++, j++) {
      if (j == r.size()) j = 0;

      point pi = r[i];
      point pj = r[j];
      
      normal_x += (((pi.z) + (pj.z)) *
		   ((pj.y) - (pi.y)));
      
      normal_y += (((pi.x) + (pj.x)) *
		   ((pj.z) - (pi.z)));
      
      normal_z += (((pi.y) + (pj.y)) *
		   ((pj.x) - (pi.x)));
    }

    return point(normal_x, normal_y, normal_z).normalize();
  }

  template<typename Ring>
  point clockwise_normal(const Ring& r) {
    return ring_normal(r);
  }

  template<typename Ring>
  point counterclockwise_normal(const Ring& r) {
    return -1*ring_normal(r);
  }
  
  template<typename Ring>
  void correct_winding_order(Ring& r, const point n) {
    double theta = angle_between(ring_normal(r), n);

    const double tol = 1.0;
    if (within_eps(theta, 0, tol)) { return; }

    if (!(within_eps(theta, 180, tol))) {
      cout << "n              = " << n << endl;
      cout << "ring_normal(r) = " << ring_normal(r) << endl;
      cout << "theta          = " << theta << endl;
      DBG_ASSERT(within_eps(theta, 180, tol));
    }

    reverse(r);

    double new_theta = angle_between(ring_normal(r), n);

    DBG_ASSERT(within_eps(new_theta, 0, tol));
  }

  template<typename Ring>
  point centroid(const Ring& r) {
    point centroid(0, 0, 0);
    for (auto p : r) {
      centroid = centroid + p;
    }

    return (1.0 / r.size()) * centroid;
  }

  template<typename Ring>
  int curve_count(const Ring& ring) {
    int count = 0;

    unsigned num_pts = ring.size();
    for (unsigned i = 0; i < num_pts; i++) {
      unsigned prev = (i + (num_pts - 1)) % num_pts;
      unsigned next = (i + 1) % num_pts;

      point prev_pt = ring[prev];
      point current_pt = ring[i];
      point next_pt = ring[next];

      point d1 = current_pt - prev_pt;
      point d2 = next_pt - current_pt;

      double angle = angle_between(d1, d2);
      cout << "Angle between = " << angle << endl;

      if (!(angle_eps(d1, d2, 0.0, 1.0) ||
	    angle_eps(d1, d2, 90.0, 3.0) ||
	    angle_eps(d1, d2, -90.0, 3.0) ||
	    angle_eps(d1, d2, 180.0, 3.0))) {
	count++;
      }
    }

    return count;
  }
  
}

#endif
