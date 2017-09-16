#include <boost/numeric/ublas/io.hpp>

#include "geometry/homogeneous_transformation.h"

namespace gca {

  homogeneous_transform
  apply(const point d, const homogeneous_transform& t) {
    return std::make_pair(t.first, t.second + to_vector(d));
  }
  
  point apply(const homogeneous_transform& t, const point p) {
    return times_3(t.first, p) + from_vector(t.second);
  }
  
  triangular_mesh apply(const homogeneous_transform& t, const triangular_mesh& m) {
    triangular_mesh rotated =
      m.apply([t](const point p)
	      { return times_3(t.first, p); });
    triangular_mesh shifted =
      rotated.apply_to_vertices([t](const point p)
				{ return p + from_vector(t.second); });
    return shifted;
  }

  boost::optional<homogeneous_transform>
  mate_planes(const plane a, const plane b, const plane c,
	      const plane ap, const plane bp, const plane cp) {

    const ublas::matrix<double> rotation =
      plane_basis_rotation(a.normal(), b.normal(), c.normal(),
			   -1*ap.normal(), -1*bp.normal(), -1*cp.normal());

    if (!within_eps(determinant(rotation), 1.0, 0.001)) {
      return boost::none;
    }

    const ublas::vector<double> displacement =
      plane_basis_displacement(rotation,
			       ap.normal(), bp.normal(), cp.normal(),
			       ap.pt(), bp.pt(), cp.pt(),
			       a.pt(), b.pt(), c.pt());

    std::pair<const ublas::matrix<double>,
	      const ublas::vector<double> > p =
      std::make_pair(rotation, displacement);
    return p;
  }

  std::vector<point> apply(const homogeneous_transform& t,
			   const std::vector<point>& pts) {
    vector<point> rpts;
    for (auto p : pts) {
      rpts.push_back(apply(t, p));
    }
    return rpts;
  }

  labeled_polygon_3 apply(const homogeneous_transform& t,
			  const labeled_polygon_3& p) {
    vector<point> rotated_verts = apply(t, p.vertices());

    vector<vector<point>> holes;
    for (auto h : p.holes()) {
      holes.push_back(apply(t, h));
    }

    polygon_3 transformed = build_clean_polygon_3(rotated_verts, holes);
    
    transformed.correct_winding_order(times_3(t.first, p.normal()));

    point rnorm = transformed.normal();
    point pnorm = p.normal();
    point rtnorm = times_3(t.first, p.normal());

    // cout << "Original normal             = " << pnorm << endl;
    // cout << "Transformed normal              = " << rnorm << endl;
    // cout << "Rotation of original normal = " << rtnorm << endl;
    
    double theta = angle_between(transformed.normal(), rtnorm);
  
    DBG_ASSERT(within_eps(theta, 0.0, 0.1));

    return transformed;
  }

}
