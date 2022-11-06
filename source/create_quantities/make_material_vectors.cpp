#include "./../../include/dynamic_fracture.h"

template <int dim>
void Dynamic_Fracture_Problem<dim>::make_material_vectors ()
{
  
  // tmp vector is used to write everything in one single vector
  // Later, this is split into 4 different vectors
  // We need 5 arguments because in material_ids are in the first column
  int no_arg = 3; // cell + mu + lame 
  Vector<double> tmp_vector(no_arg * triangulation.n_active_cells());

  // 2 vectors for the specific information on each cell
  lame_coefficient_mu_vector.reinit(triangulation.n_active_cells());
  lame_coefficient_lambda_vector.reinit(triangulation.n_active_cells());


  // Write everything from the single vector 
  // into the specific vector
  unsigned int l = 4;
  for (unsigned int i=0; i<triangulation.n_active_cells()/l; i++)
    {
      double factor = std::rand() % 500000; 
      for (unsigned int k=0; k<l; k++)
	{
	  lame_coefficient_mu_vector(l*i + k)   = factor + lame_coefficient_mu;
	  lame_coefficient_lambda_vector(l*i + k) = factor + lame_coefficient_lambda;
	}
    }
  
}

