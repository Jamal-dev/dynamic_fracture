#ifndef sneddon3D_initial_values_local_h
#define sneddon3D_initial_values_local_h
#include <vector>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
using namespace dealii;
using namespace std;


template <int dim>
  class InitialValuesSneddon3D : public Function<dim>
  {
    public:
      InitialValuesSneddon3D (const double alpha_eps) : Function<dim>(dim+1+dim) 
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

// Implementation of the initial values
  template <int dim>
  double
  InitialValuesSneddon3D<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
 
    double top = 5.0 + _alpha_eps/2.0;
    double bottom = 5.0 - _alpha_eps/2.0;
    double radius = 1.0;
    // only phase field
    if (component == dim)   
      {
	//return 1.0; 
	if ((((p(0) - 5.0) * (p(0) - 5.0) + (p(2) - 5.0) * (p(2) - 5.0)) <=  radius * radius) &&
	    ((p(1) >= bottom)  && (p(1) <= top))
	    )
	  return 0.0;
	else 
	  return 1.0;
      }
    
    return 0.0;

  }


  template <int dim>
  void
  InitialValuesSneddon3D<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesSneddon3D<dim>::value (p, comp);
  }

#endif // sneddon3D_local_h