#ifndef miehe_initial_values_local_h
#define miehe_initial_values_local_h
#include <vector>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include "./../test_cases.h"
using namespace dealii;
using namespace std;
struct PointStruct{
  double x;
  double y;
};

template <int dim>
  class InitialValuesPhaseField : public Function<dim>
  {
    public:
      InitialValuesPhaseField (const double & alpha_eps, 
			  const bool & bool_initial_crack_via_phase_field, const test_cases & test_case) : Function<dim>(dim+1+dim) 
    {
      _alpha_eps = alpha_eps;
      _bool_initial_crack_via_phase_field = bool_initial_crack_via_phase_field;
      _test_case = test_case;
      width = 10;
      height = 10;
      eps_mesh = 0.005;
      PointStruct mid_point = {width/2,height/2};
      PointStruct p1 = {0,mid_point.y-eps_mesh};
      PointStruct p2 = {0,mid_point.y+eps_mesh};
      m1 = (mid_point.y - p1.y)/(mid_point.x - p1.x);
      m2 = (mid_point.y - p2.y)/(mid_point.x - p2.x);
      b1 = mid_point.y - m1*mid_point.x;
      b2 = mid_point.y - m2*mid_point.x;
      
    }

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;

  private:
    double _alpha_eps;
    bool _bool_initial_crack_via_phase_field;
    test_cases _test_case;
    double width;
    double height;
    double eps_mesh;
    double m1,m2,b1,b2;

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
	else // _bool_initial_crack_via_phase_field == false
      {
        if (component == 2)   
          {
            
            if (_test_case.IsPMesh1())
              {
                // checking if the point lies at the crack line
                if (p(1) == m1*p(0) + b1 || p(1) == m2*p(0) + b2)
                  return 0.0;
                else
                  return 1.0;

              }            
            else
              return 1.0;
          }

      }

    // for all other components are 0
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