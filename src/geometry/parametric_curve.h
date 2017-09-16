#ifndef GCA_PARAMETRIC_CURVE_H
#define GCA_PARAMETRIC_CURVE_H

#include "utils/check.h"

namespace gca {

  class parametric_curve {
  private:
    struct concept_t {
      virtual ~concept_t() = default;
      virtual concept_t* copy() const = 0;
      virtual point value(double t) const = 0;
      virtual parametric_curve shift(point t) const = 0;
      virtual parametric_curve scale(double t) const = 0;
      virtual parametric_curve scale_xy(double t) const = 0;
      virtual parametric_curve reflect_x() const = 0;
    };

    template<typename T>
    struct model : concept_t {
      T data;
      model(T x) : data(move(x)) {}
      virtual concept_t* copy() const { return new model<T>(*this); }
      virtual point value(double t) const { return data.value(t); }
      virtual parametric_curve shift(point t) const { return data.shift(t); }
      virtual parametric_curve scale(double t) const { return data.scale(t); }
      virtual parametric_curve scale_xy(double t) const { return data.scale_xy(t); }
      virtual parametric_curve reflect_x() const { return data.reflect_x(); }
    };
  
    unique_ptr<concept_t> self;

  public:
    template<typename T>
    parametric_curve(T x) : self(new model<T>(move(x))) {}
    parametric_curve(const parametric_curve& x) : self(x.self->copy()) {}
    parametric_curve(parametric_curve&&) noexcept = default;

    parametric_curve& operator=(const parametric_curve& x)
    { parametric_curve tmp(x); *this = move(tmp); return *this; }

    parametric_curve& operator=(parametric_curve&&) noexcept = default;
    
    inline point value(double t) const
    { return self->value(t); }

    virtual parametric_curve shift(point t) const
    { return self->shift(t); }
    virtual parametric_curve scale(double t) const
    { return self->scale(t); }
    virtual parametric_curve scale_xy(double t) const
    { return self->scale_xy(t); }

    virtual parametric_curve reflect_x() const
    { return self->reflect_x(); }
    
    template<typename T>
    T& get_obj() {
      concept_t* cptr = self.get();
      parametric_curve::model<T>* mt = static_cast<parametric_curve::model<T>*>(cptr);
      return mt->data;
    }
  };

}

#endif
