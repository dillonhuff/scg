#ifndef GCA_RELATION_H
#define GCA_RELATION_H

#include <map>
#include <vector>

namespace gca {

  template<typename L, typename R>
  class relation {
  protected:
    typedef unsigned l_index;
    typedef unsigned r_index;
    
    std::vector<L> left;
    std::vector<R> right;

    std::map<l_index, std::vector<r_index> > l_to_r;
    std::map<r_index, std::vector<l_index> > r_to_l;

  public:
    relation(const std::vector<L>& ls,
	     const std::vector<R>& rs) :
      left(ls), right(rs) {
      for (unsigned i = 0; i < left.size(); i++) {
	l_to_r[i] = {};
      }

      for (unsigned i = 0; i < right.size(); i++) {
	r_to_l[i] = {};
      }
    }

    bool connected(const l_index l, const r_index r) {
      return elem(l, lefts_connected_to(r));
    }

    void insert(const l_index l, const r_index r) {
      l_to_r[l].push_back(r);
      r_to_l[r].push_back(l);
    }

    inline unsigned right_size() const { return right.size(); }

    std::vector<l_index>
    lefts_connected_to(const r_index r) const {
      auto l = r_to_l.find(r);
      assert(l != end(r_to_l));
      return l->second;
    }

    std::vector<r_index>
    rights_connected_to(const l_index l) const {
      assert(l < left_elems().size());
      auto r = l_to_r.find(l);
      assert(r != end(l_to_r));
      return r->second;
    }
    
    const std::map<l_index, std::vector<r_index> >
    left_to_right() const { return l_to_r; }

    const std::map<r_index, std::vector<l_index> >
    right_to_left() const { return r_to_l; }

    const std::vector<L>& left_elems() const { return left; }
    const std::vector<R>& right_elems() const { return right; }

    L left_elem(const unsigned i) const { return left[i]; }
    R right_elem(const unsigned i) const { return right[i]; }
    
    std::vector<unsigned> left_inds() const { return inds(left_elems()); }
    std::vector<unsigned> right_inds() const { return inds(right_elems()); }
  };
}

#endif
