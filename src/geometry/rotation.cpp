#include <boost/numeric/ublas/io.hpp>

#include "geometry/rotation.h"
#include "utils/check.h"

namespace gca {

  labeled_polygon_3 apply(const rotation& r, const labeled_polygon_3& p) {
    vector<point> pts = apply(r, p.vertices());

    vector<vector<point>> holes;
    for (auto h : p.holes()) {
      holes.push_back(apply(r, h));
    }

    polygon_3 rotated = build_clean_polygon_3(pts, holes);

    rotated.correct_winding_order(times_3(r, p.normal()));

    point rnorm = rotated.normal();
    point pnorm = p.normal();
    point rtnorm = times_3(r, p.normal());

    // cout << "Original normal             = " << pnorm << endl;
    // cout << "Rotated normal              = " << rnorm << endl;
    // cout << "Rotation of original normal = " << rtnorm << endl;
    
    double theta = angle_between(rotated.normal(), rtnorm);
  
    DBG_ASSERT(within_eps(theta, 0.0, 0.1));

    return rotated;
  }

  void test_rotation(const point from_unit, const point to_unit, const rotation& r) {
    double d = determinant(r);
    if (!(within_eps(d, 1.0, 0.001))) {
      cout << "ERROR: determinant of rotation = " << d << endl;
      cout << "from unit normal = " << from_unit << endl;
      cout << "to unit normal = " << to_unit << endl;

      double theta = angle_between(from_unit, to_unit);
      cout << "theta = " << theta << endl;
      cout << r << endl;

      point rfu = times_3(r, from_unit);
      double res_angle = angle_between(rfu, to_unit);

      cout << "r * from = " << rfu << endl;
      cout << "resulting angle = " << res_angle << endl;

      // cout << "c = " << c << endl;
      // cout << "s = " << s << endl;
      // cout << "v = " << v << endl;

      // cout << "vx = " << endl;
      // cout << vx << endl;
      
      DBG_ASSERT(false);
    }

    if (!(within_eps(angle_between(times_3(r, from_unit), to_unit), 0.0, 0.1))) {
      cout << "ERROR: Incorrect rotation " << endl;
      cout << r << endl;

      cout << "from unit normal = " << from_unit << endl;
      cout << "to unit normal = " << to_unit << endl;
      
      cout << "r*" << from_unit << " = " << times_3(r, from_unit) << " != " << to_unit << endl;

      DBG_ASSERT(false);
    }
    
  }  

  rotation rotate_unit_from_to_regular(const point from_unit, const point to_unit) {
    point v = cross(from_unit, to_unit);
    double s = v.len();
    double c = dot(from_unit, to_unit);

    boost::numeric::ublas::matrix<double> vx(3, 3);
    vx(0, 0) = 0;
    vx(0, 1) = -v.z;
    vx(0, 2) = v.y;

    vx(1, 0) = v.z;
    vx(1, 1) = 0;
    vx(1, 2) = -v.x;

    vx(2, 0) = -v.y;
    vx(2, 1) = v.x;
    vx(2, 2) = 0;

    const boost::numeric::ublas::matrix<double> id =
      boost::numeric::ublas::identity_matrix<double>(3);

    const boost::numeric::ublas::matrix<double> vxc = vx;
    
    const boost::numeric::ublas::matrix<double> vx2 = prod(vx, vx);

    const ublas::matrix<double> r =
      id + vxc + ((1.0 - c)/(s*s))*vx2;

    return r;
  }

  void check_is_unit_vector(const point v, const double tol) {
    if (!(within_eps(v.len(), 1.0, 0.00001))) {
      cout << "Error: Not a unit vector" << endl;
      cout << "v       = " << v << endl;
      cout << "v.len() = " << v.len() << endl;
      cout << "tol     = " << tol << endl;

      DBG_ASSERT(within_eps(v.len(), 1.0, tol));
    }
  }

  rotation rotate_unit_from_to(const point from_unit, const point to_unit) {

    //    DBG_ASSERT(within_eps(from_unit.len(), 1.0, 0.00001));
    //    DBG_ASSERT(within_eps(to_unit.len(),   1.0, 0.00001));

    double tol = 0.00001;

    check_is_unit_vector(from_unit, tol);
    check_is_unit_vector(to_unit, tol);
    
    double theta = angle_between(from_unit, to_unit);

    if (within_eps(theta, 0, 0.01)) {
      return boost::numeric::ublas::identity_matrix<double>(3);
    }

    if (within_eps(theta, 180, 0.01)) {
      return -1*boost::numeric::ublas::identity_matrix<double>(3);
    }

    const ublas::matrix<double> r =
      rotate_unit_from_to_regular(from_unit, to_unit);

    test_rotation(from_unit, to_unit, r);

    return r;
    
  }

  rotation rotate_from_to(const point from, const point to) {
    point from_unit = from.normalize();
    point to_unit = to.normalize();

    return rotate_unit_from_to(from_unit, to_unit);
  }

  triangular_mesh apply(const rotation& r, const triangular_mesh& m) {
    triangular_mesh rotated =
      m.apply([r](const point p)
	      { return times_3(r, p); });
    return rotated;
  }

  std::vector<point> apply(const rotation& r, const std::vector<point>& pts) {
    std::vector<point> rpts;
    for (auto p : pts) {
      rpts.push_back(times_3(r, p));
    }
    return rpts;
  }

  point apply(const rotation& r, const point p) {
    return times_3(r, p);
  }

  triangle apply(const rotation& r, const triangle tri) {
    return triangle(apply(r, tri.normal),
		    apply(r, tri.v1),
		    apply(r, tri.v2),
		    apply(r, tri.v3));
  }

  polyline apply(const rotation& r, const polyline& p) {
    vector<point> applied;
    for (auto pt : p) {
      applied.push_back(times_3(r, pt));
    }

    return polyline(applied);
  }

  boost_multipoly_2
  to_boost_multipoly_2(const rotation&r, const std::vector<polygon_3>& lines) {
    boost_multipoly_2 res;
    for (auto& pl : lines) {
      res.push_back(to_boost_poly_2(apply(r, pl)));
    }
    return res;
  }

  polygon_3 apply_no_check(const rotation& r, const polygon_3& p) {
    vector<point> pts = apply(r, p.vertices());

    vector<vector<point>> holes;
    for (auto h : p.holes()) {
      holes.push_back(apply(r, h));
    }

    polygon_3 rotated = polygon_3(pts, holes, true); //, holes);

    rotated.correct_winding_order(times_3(r, p.normal()));

    point rnorm = rotated.normal();
    point pnorm = p.normal();
    point rtnorm = times_3(r, p.normal());

    // cout << "Original normal             = " << pnorm << endl;
    // cout << "Rotated normal              = " << rnorm << endl;
    // cout << "Rotation of original normal = " << rtnorm << endl;
    
    double theta = angle_between(rotated.normal(), rtnorm);
  
    DBG_ASSERT(within_eps(theta, 0.0, 0.1));

    return rotated;
  }

  std::vector<point>
  clean_for_conversion_to_polygon_3(const std::vector<point>& vertices) {
    auto outer_ring = vertices;

    outer_ring = clean_vertices(outer_ring);

    if (outer_ring.size() < 3) { return outer_ring; }

    delete_antennas_no_fail(outer_ring);

    return outer_ring;
  }

  labeled_polygon_3
  to_labeled_polygon_3(const rotation& r, const double z, const boost_poly_2& p) {
    vector<point> vertices;
    for (auto p2d : boost::geometry::exterior_ring(p)) {
      point pt(p2d.get<0>(), p2d.get<1>(), z);
      vertices.push_back(times_3(r, pt));
    }

    vector<vector<point>> holes;
    for (auto ir : boost::geometry::interior_rings(p)) {
      vector<point> hole_verts;
      for (auto p2d : ir) {
	point pt(p2d.get<0>(), p2d.get<1>(), z);
	hole_verts.push_back(times_3(r, pt));
      }
      holes.push_back(clean_vertices(hole_verts));
    }
    return build_clean_polygon_3(clean_vertices(vertices), holes);
  }

  boost::optional<polygon_3>
  to_labeled_polygon_3_maybe(const rotation& r,
			     const double z,
			     const boost_poly_2& p) {
    vector<point> vertices;
    for (auto p2d : boost::geometry::exterior_ring(p)) {
      point pt(p2d.get<0>(), p2d.get<1>(), z);
      vertices.push_back(times_3(r, pt));
    }

    vertices = clean_for_conversion_to_polygon_3(vertices);

    if (vertices.size() < 3) { return boost::none; }

    vector<vector<point>> holes;
    for (auto ir : boost::geometry::interior_rings(p)) {
      vector<point> hole_verts;
      for (auto p2d : ir) {
	point pt(p2d.get<0>(), p2d.get<1>(), z);
	hole_verts.push_back(times_3(r, pt));
      }

      hole_verts = clean_for_conversion_to_polygon_3(hole_verts);
      if (hole_verts.size() >= 3) {
	holes.push_back(hole_verts);
      }
    }

    return build_clean_polygon_3(vertices, holes);
  }

}
