#include "catch.hpp"

#include "system/parse_stl.h"
#include "utils/arena_allocator.h"

namespace gca {

  class axis_field {
    
  };

  axis_field axis_field_for_mesh(const triangular_mesh& m) {
    axis_field f;
    return f;
  }

  void vtk_debug_axis_field(const axis_field& af) {

  }

  TEST_CASE("Axis field tests") {
    arena_allocator a;
    set_system_allocator(&a);

    triangular_mesh m =
      parse_stl("test/stl-files/onshape_parts/PSU Mount - PSU Mount.stl", 0.0001);
    
    axis_field af =
      axis_field_for_mesh(m);

    //vtk_debug_axis_field(af);
  }

}
