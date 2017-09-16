#include "geometry/offset.h"
#include "geometry/plane.h"
#include "geometry/ring.h"

namespace gca {

  double signed_distance(const plane p, const point v) {
    point u = v - p.pt();
    return signed_distance_along(u, p.normal());
  }

  boost::optional<point>
  plane_intersection(const plane p, const line l) {
    double dl = signed_distance(p, l.start);
    double dr = signed_distance(p, l.end);
    if (dl*dr > 0) {
      return boost::none;
    }
    point u = (l.start - p.pt()).normalize();
    return l.start - dl*u;
  }

  point project(const plane pl, const point p) {
    point v = project_onto(pl.pt() - p, pl.normal());
    return p + v;
  }

  std::vector<point> project(const plane pl, const std::vector<point>& pts) {
    std::vector<point> ppts;

    // TODO: Should really check that pts is planar
    for (unsigned i = 0; i < pts.size(); i++) {
      point p = pts[i];
      point pp1 = pts[(i + 1) % pts.size()];
      point diff = pp1 - p;
      if (!within_eps(angle_between(diff, pl.normal()), 0.0, 1.0) &&
      	  !within_eps(angle_between(diff, pl.normal()), 180.0, 1.0)) {
	ppts.push_back(project(pl, p));
      }
    }

    // TODO: Add tolerance as parameter?
    // TODO: Actual ring unique function?

    auto r_eq = [](const point x, const point y)
      { return components_within_eps(x, y, 0.001); };

    auto last = std::unique(begin(ppts), end(ppts), r_eq);
    ppts.erase(last, end(ppts));

    if (r_eq(ppts.front(), ppts.back())) {
      ppts.pop_back();
    }

    check_simplicity(ppts);

    return ppts;
  }

  double distance_to(const plane pl, const point p) {
    return dot(pl.normal(), pl.pt() - p);
  }

  bool right_handed(const std::vector<plane>& planes) {
    DBG_ASSERT(planes.size() == 3);

    double d = cross(planes[0].normal(), planes[1].normal()).dot(planes[2].normal());
    return d > 0.0;
  }

  std::vector<plane> set_right_handed(const std::vector<plane>& basis) {
    if (right_handed(basis)) { return basis; }
    vector<plane> rh_basis = basis;
    plane tmp = rh_basis[0];
    rh_basis[0] = rh_basis[1];
    rh_basis[1] = tmp;
    DBG_ASSERT(right_handed(rh_basis));
    return rh_basis;
  }

}
