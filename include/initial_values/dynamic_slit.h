#ifndef dynamic_slit_initial_values_local_h
#define dynamic_slit_initial_values_local_h
#include <vector>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
using namespace dealii;
using namespace std;

// this is the initial condition for the dynamic slit problem
template <int dim>
  class InitialValuesDynamicSlit : public Function<dim>
  {
    public:
      InitialValuesDynamicSlit () 
      : Function<dim>(dim+1+dim),_no_components(dim+1+dim) 
      {
      } 
   

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;

  private:
    double _no_components;

  };


  template <int dim>
  double
  InitialValuesDynamicSlit<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
   // 2 is the phase field
   if (component == 2)
    {
      return 1.0;
    }

  return 0.0;


  }


  template <int dim>
  void
  InitialValuesDynamicSlit<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesDynamicSlit<dim>::value (p, comp);
  }

#endif // miehe_initial_values_local_h