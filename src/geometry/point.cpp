#include <cassert>
#include <cmath>
#include <iostream>

#include "geometry/point.h"
#include "utils/check.h"

using namespace std;

namespace gca {

  double safe_acos(double v) {
    if (within_eps(v, -1)) { return M_PI; }
    if (within_eps(v, 1)) { return 0.0; }
    DBG_ASSERT(-1 <= v && v <= 1);
    return acos(v);
  }
  
  bool within_eps(const point& l, const point& r, double eps) {
    double xd = l.x - r.x;
    double yd = l.y - r.y;
    double zd = l.z - r.z;
    double diff = sqrt(xd*xd + yd*yd + zd*zd);
    return diff <= eps;
  }

  bool within_eps(double l, double r, double eps) {
    double diff = abs(l - r);
    return diff <= eps;
  }

  bool angle_eps(const point l, const point r, const double val, const double eps) {
    return within_eps(angle_between(l, r), val, eps);
  }

  point point::normalize() const {
    double l = len();

    //    DBG_ASSERT(!within_eps(l, 0.0));

    if (!within_eps(l, 0.0)) {
      return point(x / l, y / l, z / l);
    }

    return point(0.0, 0.0, 0.0);
  }

  
  point point::rotate_z(double degrees) const {
    double theta_rad = (M_PI/180)*degrees;
    double new_x = cos(theta_rad)*x - sin(theta_rad)*y;
    double new_y = sin(theta_rad)*x + cos(theta_rad)*y;
    return point(new_x, new_y, z);
  }

  void point::print(ostream& s) const {
    s << "(" << x << ", " << y << ", " << z << ")";
  }

  double point::len() const {
    return sqrt(x*x + y*y + z*z);
  }

  point operator*(double a, const point& p) {
    return point(a*p.x, a*p.y, a*p.z);
  }

  double angle_between(point u, point v) {
    if (within_eps(u, v)) { return 0.0; }
    double l = u.len() * v.len();
    double d = u.dot(v);
    if (within_eps(d, 0)) { return 90.0; }
    double m = (u.dot(v)) / l;
    double rads = safe_acos(m);
    return (180.0/M_PI)*rads;
  }
  
  ostream& operator<<(ostream& s, const point& p) {
    p.print(s);
    return s;
  }

  ostream& operator<<(ostream& s, const vector<point>& p) {
    for (auto pt : p) { pt.print(s); s << " "; }
    return s;
  }
  
  point extend_back(point start, point end, double l) {
    point se = end - start;
    point sp = start - ((l/se.len())*se);
    return sp;
  }

  double dot(point u, point v) {
    return u.dot(v);
  }

  point cross(point u, point v) {
    double x = u.y * v.z - v.y * u.z;
    double y = u.z * v.x - v.z * u.x;
    double z = u.x * v.y - v.x * u.y;
    return point(x, y, z);
  }

  point project_onto(point p, point proj_d) {
    point proj_dir = proj_d.normalize();
    return (p.dot(proj_dir))*proj_dir;
  }

  double signed_distance_along(const point p, const point proj_dir) {
    point dir = proj_dir.normalize();
    double len = p.dot(dir);
    return len;
  }

  double greater_than_diameter(const point normal,
			       const std::vector<point>& centroids) {
    vector<point> face_projections(centroids.size());
    transform(begin(centroids), end(centroids),
	      begin(face_projections),
	      [normal](const point cent) {
		return project_onto(cent, normal);
	      });
    auto max_e = max_element(begin(face_projections), end(face_projections),
			     [](const point l, const point r)
			     { return l.len() < r.len(); });
    double ray_len = 2*(*max_e).len();
    return ray_len;
  }

  double diameter(const point normal,
		  const std::vector<point>& pts) {
    vector<double> face_projections(pts.size());
    transform(begin(pts), end(pts),
	      begin(face_projections),
	      [normal](const point cent) {
		return signed_distance_along(cent, normal);
	      });
    auto max_e = max_element(begin(face_projections), end(face_projections));
    auto min_e = min_element(begin(face_projections), end(face_projections));
    return abs(*max_e - *min_e);
  }

  point max_along(const std::vector<point>& pts, const point normal) {
    DBG_ASSERT(pts.size() > 0);
    auto max_e = max_element(begin(pts), end(pts),
			     [normal](const point l, const point r)
			     { return signed_distance_along(l, normal) <
			       signed_distance_along(r, normal); });
    return *max_e;
  }

  point min_along(const std::vector<point>& pts, const point normal) {
    DBG_ASSERT(pts.size() > 0);
    auto min_e = min_element(begin(pts), end(pts),
			     [normal](const point l, const point r)
			     { return signed_distance_along(l, normal) <
			       signed_distance_along(r, normal); });
    return *min_e;
  }
  
  double max_distance_along(const std::vector<point>& pts, const point normal) {
    DBG_ASSERT(pts.size() > 0);
    vector<double> face_projections(pts.size());
    transform(begin(pts), end(pts),
	      begin(face_projections),
	      [normal](const point cent) {
		return signed_distance_along(cent, normal);
	      });
    auto max_e = max_element(begin(face_projections), end(face_projections));
    return *max_e;
  }

  double min_distance_along(const std::vector<point>& pts, const point normal) {
    vector<double> face_projections(pts.size());
    transform(begin(pts), end(pts),
	      begin(face_projections),
	      [normal](const point cent) {
		return signed_distance_along(cent, normal);
	      });
    auto min_e = min_element(begin(face_projections), end(face_projections));
    return *min_e;
  }

  std::vector<point> shift(const point s, const std::vector<point>& pts) {
    vector<point> res;
    for (auto p : pts) {
      res.push_back(p + s);
    }
    return res;
  }
  

  bool components_within_eps(const point l, const point r, const double tol) {
    return within_eps(l.x, r.x, tol) &&
      within_eps(l.y, r.y, tol) &&
      within_eps(l.z, r.z, tol);
  }

  bool no_duplicate_points(const std::vector<point>& pts, const double tol) {
    int num_duplicates = 0;
    for (unsigned i = 0; i < pts.size(); i++) {
      point p = pts[i];

      for (unsigned j = 0; j < pts.size(); j++) {
	if (i != j) {
	  point q = pts[j];

	  if (components_within_eps(p, q, tol)) {
	    for (unsigned k = 0; k < pts.size(); k++) {
	      cout << "pts[" << k << "] = " << pts[k] << endl;
	    }

	    cout << "DONE WITH POINTS" << endl;
	    
	    cout << "# of points = " << pts.size() << endl;
	    cout << "duplicate points at i = " << i << " , j = " << j << endl;
	    cout << "pts[i] = " << p << endl;
	    cout << "pts[j] = " << q << endl;
	    num_duplicates++;
	  }
	}
      }
    }

    if (num_duplicates > 0) {
      cout << "# of duplicate points = " << num_duplicates << endl;
    }

    return num_duplicates == 0;
  }

  
  void clean_colinear_edges(std::vector<point>& pts,
			    const double tol) {
    if (pts.size() < 3) { return; }

    bool found_colinear_edge = true;

    while (found_colinear_edge) {

      found_colinear_edge = false;

      for (unsigned i = 0; i < pts.size(); i++) {
	unsigned i1 = (i + 1) % pts.size();
	unsigned i2 = (i + 2) % pts.size();

	point p = pts[i];
	point p1 = pts[i1];
	point p2 = pts[i2];

	point raw_dir1 = p1 - p;
	point raw_dir2 = p2 - p1;

	// if (within_eps(raw_dir1.len(), 0.0)) {
	//   pts.erase(begin(pts) + i1);
	//   found_colinear_edge = true;
	//   break;
	// }

	// NOTE: This is the site of the normalization bug
	point dir1 = raw_dir1.normalize();
	point dir2 = raw_dir2.normalize();

	// point dir1 = (p1 - p).normalize();
	// point dir2 = (p2 - p1).normalize();

	if (angle_eps(dir1, dir2, 0.0, tol) ||
	    angle_eps(dir1, dir2, 180.0, tol)) {
	  pts.erase(begin(pts) + i1);
	  found_colinear_edge = true;
	  break;
	}
      }
    }
  }

  
  std::vector<point> clean_vertices(const std::vector<point>& pts,
				    const double distance_tol,
				    const double linearity_tol) {
    if (pts.size() < 3) { return pts; }

    vector<point> rpts = pts;

    bool found_duplicate = true;
    while (found_duplicate) {

      found_duplicate = false;
      
      for (unsigned i = 0; i < rpts.size(); i++) {

	point p = rpts[i];
	unsigned i1 = (i + 1) % rpts.size();
	point pp1 = rpts[i1];

	if (components_within_eps(p, pp1, distance_tol)) {
	  found_duplicate = true;
	  rpts.erase(begin(rpts) + i1);
	  break;

	}
      }

    }

    clean_colinear_edges(rpts, linearity_tol);

    return rpts;
  }

  std::vector<point> clean_vertices_within_eps(const std::vector<point>& pts,
					       const double distance_tol,
					       const double linearity_tol) {
    if (pts.size() < 3) { return pts; }

    vector<point> rpts = pts;

    bool found_duplicate = true;
    while (found_duplicate) {

      found_duplicate = false;
      
      for (unsigned i = 0; i < rpts.size(); i++) {

	point p = rpts[i];
	unsigned i1 = (i + 1) % rpts.size();
	point pp1 = rpts[i1];

	if (within_eps(p, pp1, distance_tol)) {
	  found_duplicate = true;
	  rpts.erase(begin(rpts) + i1);
	  break;

	}
      }

    }

    clean_colinear_edges(rpts, linearity_tol);

    return rpts;
  }
  
  std::vector<point> clean_vertices(const std::vector<point>& pts) {
    return clean_vertices(pts, 0.001, 0.0000001);
  }

  // TODO: Eventually move to ring
  bool has_antenna(const std::vector<point>& ring) {

    for (unsigned i = 0; i < ring.size(); i++) {
      point p = ring[i];
      point q = ring[(i + 1) % ring.size()];
      point s = ring[(i + 2) % ring.size()];
      point slope1 = q - p;
      point slope2 = s - q;

      point sum = slope1 + slope2;

      // TODO: Make this magic number tolerance an argument
      if (sum.len() < 0.001) {
	cout << "Points " << p << " " << q << " " << s << endl;
	cout << "at index " << i << " form an antenna";
	return true;
      }
    }

    return false;
  }

  void delete_antennas(std::vector<point>& pts) {
    DBG_ASSERT(pts.size() >= 3);

    bool deleted_one = true;

    while (deleted_one) {
      DBG_ASSERT(pts.size() >= 3);

      deleted_one = false;
      for (unsigned i = 0; i < pts.size(); i++) {
	point p = pts[i];
	unsigned i2 = (i + 1) % pts.size();
	unsigned i3 = (i + 2) % pts.size();
	point q = pts[i2];
	point s = pts[i3];
	point slope1 = q - p;
	point slope2 = s - q;

	point sum = slope1 + slope2;

	// TODO: Make this magic number tolerance an argument
	// TODO: Add some sort of relative magnitude test
	if (sum.len() < 0.001) {
	  cout << "Points " << p << " " << q << " " << s << endl;
	  cout << "at index " << i << " form an antenna with length = " << sum.len() << endl;;
	  cout << "in ring of size " << pts.size() << endl;
	  if (i3 > i2) {
	    pts.erase(begin(pts) + i2);
	    pts.erase(begin(pts) + i2);
	  } else {
	    pts.erase(begin(pts) + i2);
	    pts.erase(begin(pts) + i3);
	  }
	  deleted_one = true;
	}
      }
      
    }
  }

  void delete_antennas_no_fail(std::vector<point>& pts) {
    DBG_ASSERT(pts.size() >= 3);

    bool deleted_one = true;

    while (deleted_one) {
      if (pts.size() < 3) {
	return;
      }

      deleted_one = false;
      for (unsigned i = 0; i < pts.size(); i++) {
	point p = pts[i];
	unsigned i2 = (i + 1) % pts.size();
	unsigned i3 = (i + 2) % pts.size();
	point q = pts[i2];
	point s = pts[i3];
	point slope1 = q - p;
	point slope2 = s - q;

	point sum = slope1 + slope2;

	// TODO: Make this magic number tolerance an argument
	// TODO: Add some sort of relative magnitude test
	if (sum.len() < 0.001) {
	  cout << "Points " << p << " " << q << " " << s << endl;
	  cout << "at index " << i << " form an antenna with length = " << sum.len() << endl;;
	  cout << "in ring of size " << pts.size() << endl;
	  if (i3 > i2) {
	    pts.erase(begin(pts) + i2);
	    pts.erase(begin(pts) + i2);
	  } else {
	    pts.erase(begin(pts) + i2);
	    pts.erase(begin(pts) + i3);
	  }
	  deleted_one = true;
	}
      }
      
    }
  }
  
}
