#ifndef GCA_GEOMETRY_ARC_H
#define GCA_GEOMETRY_ARC_H

#include <cassert>
#include <cmath>

#include "geometry/direction.h"
#include "geometry/point.h"
#include "utils/check.h"

namespace gca {

  struct arc {
    point start, end, center;
    double radius;
    direction dir;

    arc(point startp, point endp, point start_offsetp, direction dirp) :
      start(startp), end(endp), center(start + start_offsetp),
      radius((start - center).len()), dir(dirp) {
      assert(within_eps(start.z, end.z));
      assert(within_eps(start.z, center.z));
      sanity_check();
    }

    void sanity_check() {
      double sdiff = (center - start).len();
      double ediff =  (center - end).len();
      double tolerance = 0.0005;
      if (!(within_eps(sdiff, ediff, tolerance))) {
	cout << "Error: !(within_eps((center() - start).len(), (center() - end).len()))" << endl;
	cout << "(center() - start).len() = " << sdiff << endl;
	cout << "(center() - end).len()   = " << ediff << endl;
	cout << "Difference               = " << sdiff - ediff << endl;
	cout << "Tolerance                = " << tolerance << endl;
	cout << "center                   = " << center << endl;
	cout << "start                    = " << start << endl;
	cout << "end                      = " << end << endl;
	cout << "start_offset             = " << start_offset() << endl;
	assert(false);
      }
    }

    point value(double t) const;

    inline point start_offset() const { return center - start; }

    arc shift(point s) const
    { return arc(start + s, end + s, start_offset(), dir); }

    arc reflect_x() const {
      DBG_ASSERT(false);
    }
    
    arc scale(double s) const {
      return arc(s*start, s*end, s*start_offset(), dir);
    }

    arc scale_xy(double t) const {
      point se = start;
      se.x = t*se.x;
      se.y = t*se.y;
      point ee = end;
      ee.x = t*ee.x;
      ee.y = t*ee.y;
      point so = start_offset();
      so.x = t*so.z;
      so.y = t*so.y;
      return arc(se, ee, so, dir);
    }
  };
}

#endif
