#ifndef pressurized_initial_values_local_h
#define pressurized_initial_values_local_h
#include <vector>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
using namespace dealii;
using namespace std;
template <int dim>
  class InitialValuesPressurized : public Function<dim>
  {
    public:
      InitialValuesPressurized (const double alpha_eps);

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;

  private:
    double _alpha_eps;

  };

//Implementation

template <int dim>
InitialValuesPressurized<dim>::InitialValuesPressurized (const double alpha_eps) 
: Function<dim>(dim+dim+1),
_alpha_eps(alpha_eps)
{
} 


template <int dim>
  double
  InitialValuesPressurized<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
    // only phase field
    double top    = 2.0 +    _alpha_eps;
    double bottom = 2.0 -    _alpha_eps;

    if (component == 2)   
      {
	//return 1.0; 
	if (((p(0) >= 1.5) && (p(0) <= 2.5)) &&
	    ((p(1) >= bottom) && (p(1) <= top))
	    )
	  {
	    return 0.0; 

	  }
	else 
	  return 1.0;
      }
    
    return 0.0;
  }


  template <int dim>
  void
  InitialValuesPressurized<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesPressurized<dim>::value (p, comp);
  }


#endif //pressurized_local_h