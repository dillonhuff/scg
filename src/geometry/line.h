#ifndef GCA_LINE_H
#define GCA_LINE_H

#include "geometry/point.h"

namespace gca {

  // TODO: This really belongs in utils
  template<typename T>
  struct maybe {
    bool just;
    T t;
    maybe() : just(false), t() {}
    maybe(T tp) : just(true), t(tp) {}
  };

  struct line {
    point start, end;
    line(point s, point e) : start(s), end(e) {}

    point value(double t) const
    { return (1.0 - t)*start + t*end; }

    line shift(point t) const
    { return line(start + t, end + t); }

    line scale(double t) const
    { return line(t*start, t*end); }

    line reflect_x() const {
      point reflected_start(-1*start.x, start.y, start.z);
      point reflected_end(-1*end.x, end.y, end.z);
      return line(reflected_start, reflected_end);
    }
    
    line scale_xy(double t) const {
      point se = start;
      se.x = t*se.x;
      se.y = t*se.y;
      point ee = end;
      ee.x = t*ee.x;
      ee.y = t*ee.y;
      return line(se, ee);
    }

    
  };

  bool same_line(const line l, const line r, double tolerance=0.0000001);  
  int count_in(const line l, const vector<line> ls);  
  bool adj_segment(const line l, const line r);  
  ostream& operator<<(ostream& out, line l);
  maybe<point> trim_or_extend(line prev, line next);
  point trim_or_extend_unsafe(line prev, line next);
  maybe<point> line_intersection_2d(line prev, line next);
  maybe<point> segment_intersection_2d(line p, line l);
  vector<line> make_lines(const vector<point>& pts);
}

#endif
