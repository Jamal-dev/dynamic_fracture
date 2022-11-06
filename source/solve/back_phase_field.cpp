#include "./../../include/dynamic_fracture.h"

template <int dim>
void
Dynamic_Fracture_Problem<dim>::project_back_phase_field ()  
{
  for (unsigned int i=0; i<solution.block(1).size(); ++i)
    if (solution.block(1)(i) < 0)
      solution.block(1)(i) = 0;
    else if (solution.block(1)(i) > 1)
      solution.block(1)(i) = 1;
}
