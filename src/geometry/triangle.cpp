#include <cassert>

#include "geometry/line.h"
#include "geometry/triangle.h"
#include "utils/algorithm.h"

namespace gca {

  double min_z(const vector<triangle>& triangles) {
    auto t = *min_element(begin(triangles), end(triangles),
			  [](const triangle l, const triangle r)
			  { return l.v1.z < r.v1.z; });
    return t.v1.z;
  }

  bool is_upward_facing(const triangle& t, double tolerance) {
    return (t.normal.normalize()).z > tolerance;
  }

  bool same_orientation(const triangle& x, const triangle& y, double tolerance) {
    return within_eps((x.normal - y.normal).len(), 0.0, tolerance);
  }

  vector<point> collect_polygon(vector<line>& lines) {
    DBG_ASSERT(lines.size() > 0);

    vector<point> points;
    vector<line> to_remove;
    points.push_back(lines.front().start);
    points.push_back(lines.front().end);
    lines.erase(lines.begin());
    unsigned i = 0;
    while (lines.size() > 0 && i < lines.size()) {
      if (within_eps(lines[i].start, points.back())) {
	if (within_eps(lines[i].end, points.front())) {
	  lines.erase(lines.begin() + i);

	  //DBG_ASSERT(no_duplicate_points(points, 0.001));
	  
	  return points;
	}
	points.push_back(lines[i].end);
	lines.erase(lines.begin() + i);
	i = 0;
      } else if (within_eps(lines[i].end, points.back())) {
	if (within_eps(lines[i].start, points.front())) {
	  lines.erase(lines.begin() + i);

	  //DBG_ASSERT(no_duplicate_points(points, 0.001));
	  
	  return points;
	}
	points.push_back(lines[i].start);
	lines.erase(lines.begin() + i);
	i = 0;
      } else {
	i++;
      }
    }

    DBG_ASSERT(no_duplicate_points(points, 0.001));
    
    return points;
  }

  vector<oriented_polygon>
  unordered_segments_to_polygons(point normal,
				 vector<line>& lines) {
    vector<oriented_polygon> ps;
    auto ls = lines;
    while (ls.size() > 0) {
      vector<point> vertices = collect_polygon(ls);
      ps.push_back(oriented_polygon(normal, vertices));
    }
    return ps;
  }

  vector<oriented_polygon> mesh_bounds(const vector<triangle>& triangles) {
    auto tris = triangles;
    vector<oriented_polygon> ps;
    if (tris.size() == 0) {
      return ps;
    }
    point normal = triangles.front().normal;
    vector<line> tri_lines;
    for (auto t : triangles) {
      tri_lines.push_back(line(t.v1, t.v2));
      tri_lines.push_back(line(t.v2, t.v3));
      tri_lines.push_back(line(t.v3, t.v1));
    }
    
    vector<line> no_dups;
    for (auto l : tri_lines) {
      if (count_in(l, tri_lines) == 1) {
	no_dups.push_back(l);
      }
    }

    return unordered_segments_to_polygons(normal, no_dups);
  }

  ostream& operator<<(ostream& out, const triangle& t) {
    cout << "---- TRIANGLE ----" << endl;
    cout << t.normal << endl;
    cout << t.v1 << endl;
    cout << t.v2 << endl;
    cout << t.v3 << endl;
    return out;
  }

  bool in_projection(const triangle t, const point p) {
    vector<point> vs{t.v1, t.v2, t.v3, t.v1};
    return contains(oriented_polygon(t.normal, vs), p);
  }

  double z_at(const triangle t, double x, double y) {

    point a = t.v1;
    point b = t.v2;
    point c = t.v3;
    double y_c_numerator = (b.x - a.x)*(c.z - a.z) - (c.x - a.x)*(b.z - a.z);
    double x_c_numerator = (b.y - a.y)*(c.z - a.z) - (c.y - a.y)*(b.z - a.z);
    double denom = (b.x - a.x)*(c.y - a.y) - (c.x - a.x)*(b.y - a.y);
    double y_c = y_c_numerator / denom;
    double x_c = x_c_numerator / denom;
    return a.z + y_c*(y - a.y) - x_c*(x - a.x);
  }

  bool below(const triangle t, const point p) {
    return p.z < z_at(t, p.x, p.y);
  }

  bool ray_intersects_triangle(point p, point d,
			       point v0, point v1, point v2) {
    point e1, e2, s, q;
    double a, f, u, v;
    e1 = v1 - v0;
    e2 = v2 - v0;

    point h = cross(d, e2);
    a = dot(e1, h);

    if (a > -0.00001 && a < 0.00001) {
      return(false);
    }

    f = 1/a;
    s = p - v0;
    u = f * (dot(s, h));

    if (u < 0.0 || u > 1.0)
      return(false);

    q = cross(s, e1);
    v = f * dot(d, q);

    if (v < 0.0 || u + v > 1.0) {
      return false;
    }

    // at this stage we can compute t to find out where
    // the intersection point is on the line
    double t = f * dot(e2, q);

    if (t > 0.00001) {// ray intersection
      return true;
    } else {// this means that there is a line intersection
      // but not a ray intersection
      return false;
    }

  }  

  bool intersects(const triangle t, const line l) {
    return ray_intersects_triangle(l.start, l.end - l.start,
    				   t.v1, t.v2, t.v3) &&
      ray_intersects_triangle(l.end, l.start - l.end,
			      t.v1, t.v2, t.v3);
  }

  void select_visible_triangles(vector<triangle>& triangles) {
    delete_if(triangles,
	      [](const triangle t)
	      { return !is_upward_facing(t, 0.01); });
  }

  bool intersects_triangles(line l, const vector<triangle>& triangles) {
    for (auto t : triangles) {
      if (intersects(t, l)) { return true; }
    }
    return false;
  }

  double distance_along(point normal, const triangle t) {
    point p = t.v1;
    point dir = normal.normalize();
    return ((p.dot(dir))*dir).len();
  }

  std::vector<triangle> square_triangles(const point n,
					 const std::vector<point>& pts) {
    DBG_ASSERT(pts.size() == 4);
    vector<triangle> tris;
    tris.push_back(triangle(n, pts[0], pts[1], pts[2]));
    tris.push_back(triangle(n, pts[2], pts[3], pts[0]));
    return tris;
  }

  std::vector<triangle>
  triangulate_box_pts(const std::vector<point>& p) {
    DBG_ASSERT(p.size() == 8);

    vector<point> f0{p[0], p[1], p[5], p[4]};
    vector<point> f1{p[2], p[3], p[1], p[0]};
    vector<point> f2{p[6], p[7], p[3], p[2]};
    vector<point> f3{p[4], p[5], p[7], p[6]};
    vector<point> f4{p[1], p[3], p[7], p[5]};
    vector<point> f5{p[2], p[0], p[4], p[6]};

    
    std::vector<triangle> tris;
    concat(tris, square_triangles(point(0, -1, 0), f0));
    concat(tris, square_triangles(point(-1, 0, 0), f1));
    concat(tris, square_triangles(point(0, 1, 0), f2));
    concat(tris, square_triangles(point(1, 0, 0), f3));
    concat(tris, square_triangles(point(0, 0, 1), f4));
    concat(tris, square_triangles(point(0, 0, -1), f5));
    
    DBG_ASSERT(tris.size() == 12);

    return tris;
  }

  std::vector<triangle> box_triangles(box b) {
    point n(1, 0, 0);
    point p0(b.x_min, b.y_min, b.z_min);
    point p1(b.x_min, b.y_min, b.z_max);
    point p2(b.x_min, b.y_max, b.z_min);
    point p3(b.x_min, b.y_max, b.z_max);
    point p4(b.x_max, b.y_min, b.z_min);
    point p5(b.x_max, b.y_min, b.z_max);
    point p6(b.x_max, b.y_max, b.z_min);
    point p7(b.x_max, b.y_max, b.z_max);

    return triangulate_box_pts({p0, p1, p2, p3, p4, p5, p6, p7});
  }

  unsigned find_ear_index(const std::vector<point>& pts) {
    DBG_ASSERT(pts.size() > 2);
    for (unsigned i = 1, j = 0; i < pts.size(); i++, j++) {
      unsigned ip1 = (i + 1) % pts.size();
      point l = pts[ip1];
      l.z = 0.0;
      point r = pts[j];
      r.z = 0.0;
      point mid = (1/2.0) * (l + r);
      if (contains(oriented_polygon(point(0, 0, 1), pts), mid)) {
	return i;
      }
    }
    cout << "Error: no ears left in " << endl;
    DBG_ASSERT(false);
  }

  triangle clip_ear(const unsigned ear_ind,
		    std::vector<point>& pts) {
    unsigned l = ear_ind == 0 ? pts.size() : ear_ind - 1;
    unsigned r = (ear_ind + 1) % pts.size();
    triangle t = triangle(point(0, 0, 1), pts[l], pts[ear_ind], pts[r]);
    pts.erase(begin(pts) + ear_ind);
    return t;
  }

  // TODO: Add more complex ear finding?
  std::vector<triangle> triangulate_polygon(const oriented_polygon& p) {
    DBG_ASSERT(p.vertices().size() > 2);
    vector<point> pts = p.vertices();
    vector<triangle> tris;
    while (pts.size() > 3) {
      unsigned next_ear = find_ear_index(pts);
      tris.push_back(clip_ear(next_ear, pts));
    }

    DBG_ASSERT(pts.size() == 3);
    tris.push_back(triangle(point(0, 0, 1),
			    pts[0],
			    pts[1],
			    pts[2]));

    return tris;
  }

  double max_z(const triangle t) {
    vector<double> zs = {t.v1.z, t.v2.z, t.v3.z};
    return max_e(zs);
  }

  point normal(const triangle t) {
    return cross(t.v2 - t.v1,
		 t.v3 - t.v1).normalize();
  }
}
