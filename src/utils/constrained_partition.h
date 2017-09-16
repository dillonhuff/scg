#ifndef GCA_CONSTRAINED_PARTITION_H
#define GCA_CONSTRAINED_PARTITION_H

#include "utils/partition.h"

namespace gca {

  template<typename Item, typename Bucket>
  class constrained_partition {
  protected:
    relation<Item, Bucket> constraint;
    partition<Item, Bucket> part;

  public:

    constrained_partition(const relation<Item, Bucket>& rel)
      : constraint(rel), part(rel.left_elems(), rel.right_elems()) {}
    
    // inline bool is_finished() const { return false; }

    Bucket bucket(const unsigned i) const { return part.bucket(i); }
    Item item(const unsigned i) const { return part.item(i); }

    std::vector<unsigned> bucket_inds() const { return part.bucket_inds(); }
    std::vector<unsigned> item_inds() const { return part.item_inds(); }

    // const std::vector<unsigned>& bucket_inds() const { return part.bucket_inds(); }
    // const std::vector<unsigned>& item_inds() const { return part.item_inds(); }
    
    std::vector<unsigned> items_in_bucket_inds(const unsigned i) const {
      return part.items_in_bucket_inds(i);
    }

    bool bucket_can_contain_all(const unsigned bucket_ind,
				const std::vector<unsigned>& item_inds) const {
      bool can_contain_all =
      	intersection(constraint.lefts_connected_to(bucket_ind),
      		     item_inds).size() == item_inds.size();
      return can_contain_all;
    }

    inline
    const std::vector<Item>& items() const { return part.items(); }

    inline
    const std::vector<Bucket>& buckets() const { return part.buckets(); }
    
    std::vector<unsigned> non_empty_bucket_inds() const {
      std::vector<unsigned> all_inds = bucket_inds();
      delete_if(all_inds,
    		[this](const unsigned i)
    		{ return this->items_in_bucket_inds(i).size() == 0; });
      return all_inds;
    }
    
    void assign_item_to_bucket(const unsigned i_ind,
    			       const unsigned b_ind) {
      part.assign_item_to_bucket(i_ind, b_ind);
    }

    const std::vector<unsigned> allowed_bucket_inds(const unsigned i_ind) const {
      if (!(i_ind < constraint.left_elems().size())) {
	cout << "i_ind = " << i_ind << endl;
	cout << "# items = " << constraint.left_elems().size() << endl;
	assert(false);
      }
      return constraint.rights_connected_to(i_ind);
    }

    bool item_is_assignable(const unsigned i_ind) const {
      return allowed_bucket_inds(i_ind).size() > 0;
    }

    bool all_items_are_assignable() const {
      std::vector<unsigned> i_inds = item_inds();
      return all_of(begin(i_inds), end(i_inds),
		    [this](const unsigned i_ind)
		    { return this->item_is_assignable(i_ind); });
    }

  };

}

#endif
