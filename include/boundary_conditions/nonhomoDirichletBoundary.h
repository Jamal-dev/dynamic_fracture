#ifndef non_homo_boundary_condition_local_h
#define non_homo_boundary_condition_local_h
#include <vector>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
using namespace dealii;
using namespace std;

template <int dim>
class NonhomDirichletBoundaryValues : public Function<dim> 
{
  public:
  NonhomDirichletBoundaryValues (const double time,
				 std::string test_case,
				 const double alpha_eps)    
    : Function<dim>(dim+1+dim) 
    {
      _time = time;  
      _test_case = test_case;
      _alpha_eps = alpha_eps;
    }
    
  virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const;

  virtual void vector_value (const Point<dim> &p, 
			     Vector<double>   &value) const;

private:
  double _time, _alpha_eps;
  std::string _test_case;

};




#endif