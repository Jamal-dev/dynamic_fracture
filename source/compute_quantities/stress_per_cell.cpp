#include "./../../include/dynamic_fracture.h"


template <int dim>
void Dynamic_Fracture_Problem<dim>::compute_stress_per_cell ()
{

  QGauss<dim> quad_gauss(degree+2);
  FEValues<dim> fe_values (fe, quad_gauss, update_gradients | update_quadrature_points);			   
  const unsigned int  n_q_points  = quad_gauss.size();
 std::vector<std::vector<Tensor<1,dim> > > old_solution_grads (n_q_points, 
								std::vector<Tensor<1,dim> > (dim+1+dim));
  




 typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

 //unsigned int cell_counter;

  solution_stress_per_cell.reinit(triangulation.n_active_cells());

  // This is the identity matrix in two dimensions:
  const Tensor<2,dim> Identity = ALE_Transformations
    ::get_Identity<dim> ();
 

  for (unsigned int cell_counter = 0; cell!=endc; ++cell, ++cell_counter)
    {
      fe_values.reinit (cell);
      fe_values.get_function_gradients (solution, old_solution_grads);

      double norm_stress_term = 0.0;
      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  const Tensor<2,dim> grad_u = ALE_Transformations
	    ::get_grad_u<dim> (q, old_solution_grads);

	   const Tensor<2,dim> E = 0.5 * (grad_u + transpose(grad_u));
	   const double tr_E = grad_u[0][0] + grad_u[1][1];


	  Tensor<2,dim> stress_term;
	  stress_term.clear();
	  stress_term = lame_coefficient_lambda * tr_E * Identity +  2 * lame_coefficient_mu * E;

	  norm_stress_term += ALE_Transformations
	    ::get_deviator_norm<dim> (stress_term);



	}

      solution_stress_per_cell(cell_counter) = (norm_stress_term  / n_q_points);
      

    }

}
