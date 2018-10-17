#include "geometry/point.h"

#include <cassert>
#include <iostream>

using namespace gca;
using namespace std;

int main() {
  point p(1, 1, 0);
  point q(-1, -1, 0.4);
  double a = angle_between(p, q);

  cout << "angle between " << p << " and " << q << " is " << a << endl;
  
}
