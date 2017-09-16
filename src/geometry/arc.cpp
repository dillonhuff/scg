#include "geometry/arc.h"

namespace gca {

  point arc::value(double t) const {
    point sd = start - center;
    point ed = end - center;
    double theta = angle_between(sd, ed);
    double tv = dir == COUNTERCLOCKWISE ? t : (-1)*t;
    double psi = tv*theta;
    return center + sd.rotate_z(psi);
  }

}
