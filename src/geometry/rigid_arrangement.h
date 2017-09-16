#ifndef GCA_RIGID_ARRANGEMENT_H
#define GCA_RIGID_ARRANGEMENT_H

#include <map>
#include <memory>
#include <vector>

#include "geometry/homogeneous_transformation.h"
#include "geometry/surface.h"
#include "utils/check.h"

namespace gca {

  class labeled_mesh {
  protected:
    std::map<std::string, std::vector<index_t>> labeled_surfs;
    triangular_mesh m;

  public:

    labeled_mesh(const triangular_mesh& m_p)
      : m(m_p), labeled_surfs{} {}

    labeled_mesh(const labeled_mesh& m_p)
      : m(m_p.mesh()), labeled_surfs(m_p.labeled_surfs) {}
    
    surface labeled_surface(const std::string& name) const {
      auto r = labeled_surfs.find(name);
      if (r == end(labeled_surfs)) {
	cout << "Could not find " << name << endl;
	cout << "Labels:" << endl;
	for (auto n : labeled_surfs) {
	  cout << n.first << endl;
	}
	DBG_ASSERT(false);
      }

      return surface(&m, r->second);
    }

    triangular_mesh& mesh() { return m; }
    const triangular_mesh& mesh() const { return m; }

    void insert_label(const std::string& label_name,
		      const std::vector<index_t>& triangle_inds) {
      DBG_ASSERT(labeled_surfs.find(label_name) == end(labeled_surfs));
      labeled_surfs[label_name] = triangle_inds;
    }

    labeled_mesh apply(const homogeneous_transform& t) const {
      labeled_mesh m(gca::apply(t, mesh()));
      m.labeled_surfs = labeled_surfs;
      return m;
    }
    

  };

  struct arrangement_metadata {
    bool display_during_debugging;

    arrangement_metadata()
      : display_during_debugging(true) {}
  };

  class rigid_arrangement {
    std::vector<std::unique_ptr<labeled_mesh>> meshes;
    std::map<std::string, labeled_mesh*> name_index;
    std::map<std::string, arrangement_metadata> metadata_index;
    
  public:

    rigid_arrangement() {}

    rigid_arrangement(const rigid_arrangement& x) {
      for (auto n : x.mesh_names()) {
	insert(n, x.labeled_mesh(n));
	set_metadata(n, x.metadata(n));
      }
    }

    rigid_arrangement(rigid_arrangement&&) noexcept = default;

    rigid_arrangement& operator=(const rigid_arrangement& x)
    { rigid_arrangement tmp(x); *this = std::move(tmp); return *this; }

    rigid_arrangement& operator=(rigid_arrangement&&) noexcept = default;
    
    void insert(const std::string& name,
		const triangular_mesh& m) {
      meshes.push_back(unique_ptr<class labeled_mesh>(new class labeled_mesh(triangular_mesh(m))));
      std::unique_ptr<class labeled_mesh>* res = &(meshes.back());
      name_index[name] = res->get();
      metadata_index[name] = arrangement_metadata();
    }

    bool contains_mesh(const std::string& name) const
    { return name_index.find(name) != end(name_index); }

    void insert(const std::string& name,
    		const class labeled_mesh& m) {
      meshes.push_back(unique_ptr<class labeled_mesh>(new class labeled_mesh(m)));
      std::unique_ptr<class labeled_mesh>* res = &(meshes.back());
      name_index[name] = res->get();
      metadata_index[name] = arrangement_metadata();
    }
    
    void insert(const std::string& name,
		const homogeneous_transform& t,
		const triangular_mesh& m) {
      DBG_ASSERT(name_index.find(name) == end(name_index));
      insert(name, apply(t, m));
    }

    void insert(const std::string& name,
		const homogeneous_transform& t,
		const gca::labeled_mesh& m) {
      DBG_ASSERT(name_index.find(name) == end(name_index));
      insert(name, m.apply(t));
    }
    
    const std::vector<std::string> mesh_names() const {
      std::vector<std::string> names;
      for (auto p : name_index) {
	names.push_back(p.first);
      }
      return names;
    }

    const labeled_mesh& labeled_mesh(const std::string& name) const {
      auto r = name_index.find(name);
      DBG_ASSERT(r != end(name_index));
      return *(r->second);
    }

    gca::labeled_mesh& labeled_mesh(const std::string& name) {
      auto r = name_index.find(name);
      DBG_ASSERT(r != end(name_index));
      return *(r->second);
    }
    
    const triangular_mesh& mesh(const std::string& name) const {
      auto r = name_index.find(name);
      DBG_ASSERT(r != end(name_index));
      return r->second->mesh();
    }

    triangular_mesh& mesh(const std::string& name) {
      auto r = name_index.find(name);
      DBG_ASSERT(r != end(name_index));
      return r->second->mesh();
    }
    

    arrangement_metadata& metadata(const std::string& name) {
      auto r = metadata_index.find(name);
      DBG_ASSERT(r != end(metadata_index));
      return (r->second);
    }

    void set_metadata(const std::string& name, const arrangement_metadata& new_meta) {
      auto r = metadata_index.find(name);
      DBG_ASSERT(r != end(metadata_index));
      metadata_index[name] = new_meta;
    }
    
    const arrangement_metadata& metadata(const std::string& name) const {
      auto r = metadata_index.find(name);
      DBG_ASSERT(r != end(metadata_index));
      return (r->second);
    }

    surface labeled_surface(const std::string& mesh_name,
			    const std::string& surface_name) const {

      auto& s = labeled_mesh(mesh_name);

      return s.labeled_surface(surface_name);
    }

    void insert_label(const std::string& mesh_name,
		      const std::string& label_name,
		      const std::vector<index_t>& triangle_inds) {
      gca::labeled_mesh& s = labeled_mesh(mesh_name);
      s.insert_label(label_name, triangle_inds);
    }

  };

}

#endif
