#ifndef miehe_initial_values_local_h
#define miehe_initial_values_local_h
#include <vector>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
using namespace dealii;
using namespace std;


template <int dim>
  class InitialValuesPhaseField : public Function<dim>
  {
    public:
      InitialValuesPhaseField (const double alpha_eps, 
			  const bool bool_initial_crack_via_phase_field) : Function<dim>(dim+1+dim) 
    {
      _alpha_eps = alpha_eps;
_bool_initial_crack_via_phase_field = bool_initial_crack_via_phase_field;
    }

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;

  private:
    double _alpha_eps;
    bool _bool_initial_crack_via_phase_field;

  };


  template <int dim>
  double
  InitialValuesPhaseField<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
    // only phase field
	if (_bool_initial_crack_via_phase_field)
	{
    double top = 0.5 + _alpha_eps/2.0;
    double bottom = 0.5 - _alpha_eps/2.0;

    if (component == 2)   
      {
	//return 1.0; 
	if (((p(0) >= 0.5) && (p(0) <= 1.0)) &&
	    ((p(1) >= bottom) && (p(1) <= top))
	    )
	  {
	    return 0.0; 

	  }
	else 
	  return 1.0;
      }
}
	else
{
 if (component == 2)   
      {
	return 1.0;
      }

}

    
    return 0.0;


  }


  template <int dim>
  void
  InitialValuesPhaseField<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesPhaseField<dim>::value (p, comp);
  }

#endif // miehe_initial_values_local_h