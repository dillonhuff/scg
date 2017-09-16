#include <cassert>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "geometry/matrix.h"

namespace gca {

  // matrix<3, 3> rotate_onto(point a_v, point b_v) {
  //   point a = a_v.normalize();
  //   point b = b_v.normalize();
  //   point v = cross(a, b);
  //   double s = v.len();
  //   double c = a.dot(b);
  //   if (within_eps(s, 0)) {
  //     return identity<3, 3>();
  //   }
  //   double cs = ((1 - c) / (s*s));
  //   matrix<3, 3> id = identity<3, 3>();
  //   double v_coeffs[] = {0, -v.z, v.y, v.z, 0, -v.x, -v.y, v.x, 0};
  //   matrix<3, 3> vx(v_coeffs);
  //   return id + vx + cs*(vx*vx);
  // }

  // point operator*(const matrix<3, 3>& m, const point a) {
  //   double xr = m.get(0, 0) * a.x + m.get(0, 1) * a.y + m.get(0, 2) * a.z;
  //   double yr = m.get(1, 0) * a.x + m.get(1, 1) * a.y + m.get(1, 2) * a.z;
  //   double zr = m.get(2, 0) * a.x + m.get(2, 1) * a.y + m.get(2, 2) * a.z;
  //   point res(xr, yr, zr);
  //   return res;
  // }

  int determinant_sign(const ublas::permutation_matrix<std ::size_t>& pm) {
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
      if (i != pm(i))
	pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
  }

  double determinant(ublas::matrix<double>& m) {
    ublas::permutation_matrix<size_t> pm(m.size1());
    double det = 1.0;
    if( lu_factorize(m,pm) ) {
      det = 0.0;
    } else {
      for(int i = 0; i < m.size1(); i++)
	det *= m(i,i); // multiply by elements on diagonal
      det = det * determinant_sign( pm );
    }
    return det;
  }

  double determinant(const ublas::matrix<double>& m) {
    ublas::matrix<double> l = m;
    return determinant(l);
  }

  ublas::matrix<double> inverse(ublas::matrix<double>& a) {
    ublas::matrix<double> a_inv = ublas::identity_matrix<double>(a.size1());
    ublas::permutation_matrix<size_t> pm(a.size1());
    int res = lu_factorize(a, pm);
    if (!res) {
      lu_substitute(a, pm, a_inv);
    } else {
      std::cout << a << std::endl;
      std::cout << "Singular matrix!" << std::endl;
      assert(false);
    }
    return a_inv;
  }

  ublas::matrix<double> inverse(const ublas::matrix<double>& a) {
    ublas::matrix<double> b = a;
    return inverse(b);
  }

  using namespace gca;

  ublas::matrix<double>
  plane_basis_rotation(const point at, const point bt, const point ct,
		       const point apt, const point bpt, const point cpt) {
    boost::numeric::ublas::matrix<double> a(3, 3);
    a(0, 0) = at.x;
    a(1, 0) = at.y;
    a(2, 0) = at.z;

    a(0, 1) = bt.x;
    a(1, 1) = bt.y;
    a(2, 1) = bt.z;

    a(0, 2) = ct.x;
    a(1, 2) = ct.y;
    a(2, 2) = ct.z;
  
    boost::numeric::ublas::matrix<double> b(3, 3);
    b(0, 0) = apt.x;
    b(1, 0) = apt.y;
    b(2, 0) = apt.z;

    b(0, 1) = bpt.x;
    b(1, 1) = bpt.y;
    b(2, 1) = bpt.z;

    b(0, 2) = cpt.x;
    b(1, 2) = cpt.y;
    b(2, 2) = cpt.z;

    auto a_inv = inverse(a);
    return prod(b, a_inv);
  }

  boost::numeric::ublas::vector<double>
  to_vector(const point p) {
    boost::numeric::ublas::vector<double> v(3);
    v(0) = p.x;
    v(1) = p.y;
    v(2) = p.z;
    return v;
  }

  point from_vector(ublas::vector<double> v) {
    return point(v(0), v(1), v(2));
  }

  boost::numeric::ublas::vector<double>
  plane_basis_displacement(const boost::numeric::ublas::matrix<double>& r,
			   const point u1, const point u2, const point u3,
			   const point q1, const point q2, const point q3,
			   const point p1, const point p2, const point p3) {
    boost::numeric::ublas::vector<double> uv1 = to_vector(u1);
    boost::numeric::ublas::vector<double> uv2 = to_vector(u2);
    boost::numeric::ublas::vector<double> uv3 = to_vector(u3);

    boost::numeric::ublas::vector<double> qv1 = to_vector(q1);
    boost::numeric::ublas::vector<double> qv2 = to_vector(q2);
    boost::numeric::ublas::vector<double> qv3 = to_vector(q3);

    boost::numeric::ublas::vector<double> pv1 = to_vector(p1);
    boost::numeric::ublas::vector<double> pv2 = to_vector(p2);
    boost::numeric::ublas::vector<double> pv3 = to_vector(p3);

    boost::numeric::ublas::vector<double> s(3);
    s(0) = inner_prod(qv1, uv1) - inner_prod(prod(r, pv1), uv1);
    s(1) = inner_prod(qv2, uv2) - inner_prod(prod(r, pv2), uv2);
    s(2) = inner_prod(qv3, uv3) - inner_prod(prod(r, pv3), uv3);

    boost::numeric::ublas::matrix<double> u(3, 3);
    u(0, 0) = uv1(0);
    u(0, 1) = uv1(1);
    u(0, 2) = uv1(2);

    u(1, 0) = uv2(0);
    u(1, 1) = uv2(1);
    u(1, 2) = uv2(2);

    u(2, 0) = uv3(0);
    u(2, 1) = uv3(1);
    u(2, 2) = uv3(2);

    auto u_inv = inverse(u);

    return prod(u_inv, s);
  }

  point times_3(const ublas::matrix<double> m, const point p) {
    auto v = to_vector(p);
    return from_vector(prod(m, v));
  }
}
