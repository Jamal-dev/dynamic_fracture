#ifndef het3D_initial_values_local_h
#define het3D_initial_values_local_h

#include <deal.II/lac/vector.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
using namespace dealii;
using namespace std;


template <int dim>
  class InitialValuesHet3D : public Function<dim>
  {
    public:
      InitialValuesHet3D (const double alpha_eps) : Function<dim>(dim+1+dim) 
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
  InitialValuesHet3D<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
 
    double width = _alpha_eps/2.0;
    double top = 5.0 + _alpha_eps/2.0;
    double bottom = 5.0 - _alpha_eps/2.0;
    double radius = 1.0;
    // only phase field
    if (component == dim)   
      {
	if (((p(0) >= 2.6 - width/2.0) && (p(0) <= 2.6 + width/2.0))
	    && ((p(1) >= 3.8- width/2.0) && (p(1) <= 5.5 + width/2.0))
	    && (p(2) >=4 - width/2.0) && (p(2) <=4 + width/2.0)
	    )
	  return 0.0;
	else if (((p(0) >= 5.5) && (p(0) <= 7.0))
		 && ((p(1) >= 4.0 - width/2.0) && (p(1) <= 4.0 + width/2.0))
		 && (p(2) >=6 - width/2.0) && (p(2) <=6 + width/2.0)
		 )
	  return 0.0;
	else
	  return 1.0;


      }
    
    return 0.0;

  }


  template <int dim>
  void
  InitialValuesHet3D<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesHet3D<dim>::value (p, comp);
  }

#endif // het3D_local_h