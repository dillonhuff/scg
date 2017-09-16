#pragma once

#include "utils/check.h"

namespace gca {

  template<typename T, typename Eq>
  std::vector<std::vector<T> > group_by(const std::vector<T>& elems,
					Eq eq) {
    std::vector<std::vector<T> > groups;
    for (auto& e : elems) {

      bool found_group = false;
      for (auto& g : groups) {
	if (eq(e, g.front())) {
	  found_group = true;
	  g.push_back(e);
	  break;
	}
      }

      if (!found_group) {
	groups.push_back({e});
      }
    }

    return groups;
  }

}
