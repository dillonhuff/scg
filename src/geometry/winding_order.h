#pragma once

namespace gca {

  template<typename Triangle>
  bool
  share_edge(const Triangle tl,
	     const Triangle tr) {
    int num_eq = 0;
    for (unsigned i = 0; i < 3; i++) {
      for (unsigned j = 0; j < 3; j++) {
	num_eq += (get_vertex(tl, i) == get_vertex(tr, j)) ? 1 : 0;
      }
    }
    return num_eq > 1;
  }

  
  template<typename Triangle>
  bool winding_conflict(const Triangle ti, const Triangle tj) {
    for (unsigned l = 0; l < 3; l++) {
      unsigned lp1 = (l + 1) % 3;
      for (unsigned k = 0; k < 3; k++) {
	unsigned kp1 = (k + 1) % 3;

	if (get_vertex(ti, k) == get_vertex(tj, l) &&
	    get_vertex(ti, kp1) == get_vertex(tj, lp1)) {
	  return true;
	}

      }
    }
    return false;
  }

  template<typename Triangle>
  Triangle flip_winding_order(const Triangle& t) {
    Triangle f;
    set_vertex(f, 0, get_vertex(t, 1));
    set_vertex(f, 1, get_vertex(t, 0));
    set_vertex(f, 2, get_vertex(t, 2));
    
    return f;
  }

  template<typename Triangle>
  std::vector<Triangle>
  flip_winding_orders(const std::vector<Triangle>& vertex_triangles) {
    vector<Triangle> tris;
    for (auto t : vertex_triangles) {
      tris.push_back(flip_winding_order(t));
    }
    return tris;
  }

  template<typename Triangle>
  Triangle
  correct_orientation(const Triangle to_correct,
		      const std::vector<Triangle>& others) {
    auto ti = to_correct;
    for (auto tj : others) {
      if (winding_conflict(ti, tj)) {
	Triangle corrected = flip_winding_order(ti);

	return corrected;
      }

    }

    return ti;
  }

  template<typename Triangle>
  int
  num_winding_order_errors(const std::vector<Triangle>& triangles) {
    int num_errs = 0;
    for (unsigned i = 0; i < triangles.size(); i++) {
      for (unsigned j = i; j < triangles.size(); j++) {
	if (i != j) {
	  auto ti = triangles[i];
	  auto tj = triangles[j];

	  if (winding_conflict(ti, tj)) {
	    num_errs++;
	  }
	  
	}
      }
    }
    return num_errs;
  }

  template<typename Triangle>
  std::vector<Triangle>
  fix_wind_errors(const std::vector<Triangle>& triangles) {
    vector<Triangle> tris;
    vector<unsigned> remaining_inds = inds(triangles);

    unsigned num_added = 0;
    
    while (remaining_inds.size() > 0) {

      for (auto ind : remaining_inds) {
	Triangle next_t = triangles[ind];
	vector<Triangle> sub_tris =
	  select(tris, [next_t](const Triangle t)
		 { return share_edge(next_t, t); });

	if (sub_tris.size() > 0) {
	  Triangle corrected = correct_orientation(next_t, sub_tris);
	  tris.push_back(corrected);

	  remove(ind, remaining_inds);
	  num_added++;
	  break;
	}

	if (num_added == 0) {
	  tris.push_back(next_t);

	  remove(ind, remaining_inds);
	  num_added++;
	  break;
	}
      }

    }


    return tris;
  }
  
  // TODO: Remove use of utils/algorithm
  template<typename Triangle>
  std::vector<Triangle>
  fix_winding_order_errors(const std::vector<Triangle>& triangles) {

    vector<Triangle> tris = fix_wind_errors(triangles);

    auto ccs =
      connected_components_by(tris, [](const Triangle l, const Triangle r)
			      { return share_edge(l, r); });
    DBG_ASSERT(ccs.size() == 1);

    auto num_errs_after_correction = num_winding_order_errors(tris);
    if (num_errs_after_correction != 0) {
      cout << "ERROR: " << num_errs_after_correction << " winding errors still exist" << endl;
      DBG_ASSERT(num_errs_after_correction == 0);
    }

    DBG_ASSERT(num_winding_order_errors(tris) == 0);
    DBG_ASSERT(tris.size() == triangles.size());

    return tris;
  }

  
}
