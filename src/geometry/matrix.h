#ifndef GCA_MATRIX_H
#define GCA_MATRIX_H

#include <cstdlib>
#include <boost/numeric/ublas/matrix.hpp>

#include "geometry/point.h"

namespace ublas = boost::numeric::ublas;

namespace gca {

  typedef ublas::matrix<double> matrix;
  typedef ublas::vector<double> vec;

  point from_vector(ublas::vector<double> v);
  ublas::vector<double> to_vector(const point p);

  double determinant(matrix& m);
  double determinant(const matrix& m);

  matrix inverse(matrix& a);
  matrix inverse(const matrix& a);

  matrix
  plane_basis_rotation(const point at, const point bt, const point ct,
		       const point apt, const point bpt, const point cpt);

  boost::numeric::ublas::vector<double>
  plane_basis_displacement(const matrix& r,
			   const point u1, const point u2, const point u3,
			   const point q1, const point q2, const point q3,
			   const point p1, const point p2, const point p3);
  
  point times_3(const matrix m, const point p);

}

#endif
