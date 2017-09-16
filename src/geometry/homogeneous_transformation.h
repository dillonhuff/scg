#ifndef GCA_HOMOGENEOUS_TRANSFORMATION_H
#define GCA_HOMOGENEOUS_TRANSFORMATION_H

#include <boost/optional.hpp>

#include "geometry/matrix.h"
#include "geometry/plane.h"
#include "geometry/polygon_3.h"
#include "geometry/triangular_mesh.h"

namespace gca {

  typedef std::pair<const ublas::matrix<double>, const ublas::vector<double>>
    homogeneous_transform;

  triangular_mesh apply(const homogeneous_transform& t, const triangular_mesh& m);

  boost::optional<homogeneous_transform>
  mate_planes(const plane a, const plane b, const plane c,
	      const plane ap, const plane bp, const plane cp);

  point apply(const homogeneous_transform& t, const point p);

  homogeneous_transform
  apply(const point d, const homogeneous_transform& t);

    polygon_3 apply(const homogeneous_transform& r, const polygon_3& p);

}

#endif
