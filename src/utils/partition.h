#ifndef GCA_PARTITION_H
#define GCA_PARTITION_H

#include "utils/relation.h"

namespace gca {

  template<typename Item, typename Bucket>
  class partition {
  protected:
    relation<Item, Bucket> rel;

  public:

    partition(const std::vector<Item>& p_items)
      : rel{p_items, {}} {}

    partition(const std::vector<Item>& p_items,
	      const std::vector<Bucket>& p_buckets)
      : rel{p_items, p_buckets} {}
    
    inline bool is_finished() const { return false; }

    Bucket bucket(const unsigned i) const { return rel.right_elem(i); }
    Item item(const unsigned i) const { return rel.left_elem(i); }

    std::vector<unsigned> bucket_inds() const { return rel.right_inds(); }
    std::vector<unsigned> item_inds() const { return rel.left_inds(); }

    std::vector<unsigned> items_in_bucket_inds(const unsigned i) const {
      return rel.lefts_connected_to(i);
    }

    inline
    const std::vector<Item>& items() const { return rel.left_elems(); }

    inline
    const std::vector<Bucket>& buckets() const { return rel.right_elems(); }
    
    std::vector<unsigned> non_empty_bucket_inds() const {
      std::vector<unsigned> all_inds = bucket_inds();
      delete_if(all_inds,
		[this](const unsigned i)
		{ return this->items_in_bucket_inds(i).size() == 0; });
      return all_inds;
    }
    
    void assign_item_to_bucket(const unsigned i_ind,
			       const unsigned b_ind) {
      rel.insert(i_ind, b_ind);
    }

  };

}

#endif
