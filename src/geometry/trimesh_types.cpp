#include "geometry/trimesh_types.h"

namespace gca {

  bool operator==(const triangle_t l, const triangle_t r) {
    return (l.v[0] == r.v[0]) && (l.v[1] == r.v[1]) && (l.v[2] == r.v[2]);
  }

  bool operator!=(const triangle_t l, const triangle_t r) {
    return !(l == r);
  }
}
