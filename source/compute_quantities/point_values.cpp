#include "./../../include/dynamic_fracture.h"


// With help of this function, we extract 
// point values for a certain component from our
// discrete solution. We use it to gain the 
// displacements of the structure in the x- and y-directions.
template <int dim>
double Dynamic_Fracture_Problem<dim>::compute_point_value (Point<dim> p, 
					       const unsigned int component) const  
{
 
  Vector<double> tmp_vector(dim+1);
  VectorTools::point_value (dof_handler, 
			    solution, 
			    p, 
			    tmp_vector);
  
  return tmp_vector(component);
}
