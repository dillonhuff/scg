#ifndef GCA_POINT_H
#define GCA_POINT_H

#include <iostream>
#include <vector>

#include "geometry/point.h"

using namespace std;

namespace gca {

  class point {
  public:
    double x, y, z;

    point() {
      x = 0;
      y = 0;
      z = 0;
    }

  point(double xp, double yp, double zp) :
    x(xp), y(yp), z(zp) {}

    bool operator==(const point& other) const {
      return x == other.x && y == other.y && z == other.z;
    }

    bool operator!=(const point& other) const {
      return !(*this == other);
    }
    
    point operator+(const point& other) const {
      return point(x + other.x, y + other.y, z + other.z);
    }

    point operator-(const point& other) const {
      return point(x - other.x, y - other.y, z - other.z);
    }

    double len() const;

    point normalize() const;

    point rotate_z(double degrees) const;

    void print(ostream& s) const;

    double dot(point v) const {
      return x*v.x + y*v.y + z*v.z;
    }

  };

  double safe_acos(double v);
  
  point cross(point b, point c);
  double dot(point u, point v);

  point operator*(double a, const point& other);
  bool within_eps(const point& l, const point& r, double eps=0.0000001);
  bool within_eps(double l, double r, double eps=0.0000001);
  point extend_back(point start, point end, double l);
  double angle_between(point u, point v);
  std::ostream& operator<<(std::ostream& s, const point& p);
  std::ostream& operator<<(std::ostream& s, const std::vector<point>& p);
  point project_onto(point p, point proj_d);

  double signed_distance_along(const point p, const point proj_dir);

  double max_distance_along(const std::vector<point>& pts, const point proj_dir);
  double min_distance_along(const std::vector<point>& pts, const point proj_dir);

  point max_along(const std::vector<point>& pts, const point proj_dir);
  point min_along(const std::vector<point>& pts, const point proj_dir);

  std::vector<point> shift(const point s, const std::vector<point>& pts);
  double greater_than_diameter(const point normal,
			       const std::vector<point>& centroids);

  double diameter(const point normal,
		  const std::vector<point>& pts);

  bool components_within_eps(const point l, const point r, const double tol);

  bool no_duplicate_points(const std::vector<point>& pts, const double tol);

  std::vector<point> clean_vertices(const std::vector<point>& pts);

  bool angle_eps(const point l, const point r, const double val, const double eps);

  void clean_colinear_edges(std::vector<point>& pts,
			    const double tol);

  std::vector<point> clean_vertices_within_eps(const std::vector<point>& pts,
					       const double distance_tol,
					       const double linearity_tol);

  void delete_antennas(std::vector<point>& pts);

  void delete_antennas_no_fail(std::vector<point>& pts);

  bool has_antenna(const std::vector<point>& ring);

}

#endif
