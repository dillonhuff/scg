#include "geometry/point.h"
#include "geometry/rotation.h"
#include "system/parse_stl.h"

#include <cassert>
#include <iostream>


using namespace gca;
using namespace std;

void angle_between_points() {
  point p(1, 1, 0);
  point q(-1, -1, 0.4);
  double a = angle_between(p, q);

  cout << "angle between " << p << " and " << q << " is " << a << endl;
}


void rotate_onto() {
  point from(1, 0, 0);
  point to(0, 0, 1);

  const rotation r = rotate_from_to(from, to);

  cout << apply(r, from) << endl;
}

void load_mesh() {
  auto triangles = parse_stl("./test/stl-files/SlicedCone.stl").triangles;
  triangular_mesh m = make_mesh(triangles, 0.001);
  cout << "# of triangles = " << m.triangle_list().size() << endl;
}

int main() {
  angle_between_points();
  rotate_onto();
  load_mesh();
}
