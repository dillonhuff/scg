#ifndef GCA_OFFSET_H
#define GCA_OFFSET_H

#include <boost/shared_ptr.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include "geometry/polygon.h"
#include "geometry/polygon_3.h"
#include "utils/check.h"

namespace gca {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_2                    Point;
  typedef CGAL::Polygon_2<K>            Polygon_2;

  Polygon_2
  CGAL_polygon_for_oriented_polygon(const oriented_polygon& p);

  void check_simplicity(const oriented_polygon& p);
  
  std::vector<oriented_polygon> exterior_offset(const oriented_polygon& p,
						const double inc);

  std::vector<oriented_polygon> interior_offset(const oriented_polygon& p,
						const double inc);

  void check_simplicity(const std::vector<point>& pts);

  boost::optional<polygon_3>
  clean_polygon_for_offsetting_maybe(const polygon_3& poly);

  std::vector<point>
  clean_ring_for_offsetting_no_fail(const std::vector<point>& ring);

  std::vector<point>
  exterior_offset(const std::vector<point>& pts, const double tol);

  std::vector<point>
  interior_offset(const std::vector<point>& pts, const double tol);

  std::vector<std::vector<point>>
  interior_offsets(const std::vector<point>& pts,
		   const double tol);

  polygon_3 exterior_offset(const polygon_3& polys,
			    const double d);
  
  std::vector<polygon_3> exterior_offset(const std::vector<polygon_3>& polys,
					 const double d);
  std::vector<polygon_3> interior_offset(const std::vector<polygon_3>& polys,
					 const double d);

}

#endif
