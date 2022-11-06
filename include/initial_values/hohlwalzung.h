#ifndef hohlwalzung_initial_values_local_h
#define hohlwalzung_initial_values_local_h
#include <vector>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
using namespace dealii;
using namespace std;


template <int dim>
  class InitialValuesHohlwalzung : public Function<dim>
  {
    public:
      InitialValuesHohlwalzung (const double alpha_eps) : Function<dim>(dim+1+dim) 
    {
      _alpha_eps = alpha_eps;
    }

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;

  private:
    double _alpha_eps;

  };


  template <int dim>
  double
  InitialValuesHohlwalzung<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
  
    if (component == dim)
      {
	if (((p(0) >= 0.0 - _alpha_eps) && (p(0) <= 0.0 + _alpha_eps)) &&
	    ((p(1) >= (-13.0)) && (p(1) <= (-7.0))) // -13.0 => 6mm
	    )
	  return 0.0;
	else 
	  return 1.0;
      }
    

    return 0.0;
  }


  template <int dim>
  void
  InitialValuesHohlwalzung<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesHohlwalzung<dim>::value (p, comp);
  }

#endif // hohlwalzung_initial_values_local_h