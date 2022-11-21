#ifndef non_homo_boundary_condition_velocity_local_h
#define non_homo_boundary_condition_velocity_local_h
#include <vector>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include "./../test_cases.h"
using namespace dealii;
using namespace std;

template <int dim>
class NonhomDirichletBoundaryVelocity : public Function<dim> 
{
  public:
  NonhomDirichletBoundaryVelocity (const double time,
				 test_cases test_case,
				 const double alpha_eps)    
    : Function<dim>(dim+1+dim) 
    {
      _time = time;  
      _test_case = test_case;
      _alpha_eps = alpha_eps;
    }
    
  struct components {unsigned int disp_x =0 , disp_y = 1; 
                     unsigned int phase_field = 2; 
                     unsigned int vel_x = 3, vel_y=4;}comp;
  virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const;

  virtual void vector_value (const Point<dim> &p, 
			     Vector<double>   &value) const;

private:
  double _time, _alpha_eps;
  std::string _test_case;

};




#endif