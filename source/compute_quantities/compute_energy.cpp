#include "./../../include/dynamic_fracture.h"

template <int dim>
void
Dynamic_Fracture_Problem<dim>::compute_energy()
{
  // What are we computing? In Latex-style it is:
  // bulk energy = [(1+k)phi^2 + k] psi(e)
  // crack energy = \frac{G_c}{2}\int_{\Omega}\Bigl( \frac{(\varphi - 1)^2}{\eps}
  //+ \eps |\nabla \varphi|^2 \Bigr) \, dx

  const QGauss<dim> quadrature_formula(degree+2);
  const unsigned int n_q_points = quadrature_formula.size();

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_quadrature_points | update_JxW_values
                          | update_gradients);

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active(), endc = dof_handler.end();

 std::vector<Vector<double> > 
    old_solution_values (n_q_points, Vector<double>(dim+1+dim));

  std::vector<std::vector<Tensor<1,dim> > > 
    old_solution_grads (n_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));


  double local_bulk_energy = 0.0;
  double local_crack_energy = 0.0;

  for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);


        fe_values.get_function_values(solution, old_solution_values);
        fe_values.get_function_gradients(solution, old_solution_grads);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
	    const Tensor<2,dim> grad_u = ALE_Transformations 
	      ::get_grad_u<dim> (q, old_solution_grads);

            const Tensor<2,dim> E = 0.5 * (grad_u + transpose(grad_u));
            const double tr_E = dealii::trace(E);

	    const double pf = old_solution_values[q](dim);

            const double tr_e_2 = dealii::trace(E*E);

            const double psi_e = 0.5 * lame_coefficient_lambda * tr_E*tr_E + lame_coefficient_mu * tr_e_2;

            local_bulk_energy += ((1+constant_k)*pf*pf+constant_k) * psi_e * fe_values.JxW(q);

            local_crack_energy += G_c/2.0 * ((pf-1) * (pf-1)/alpha_eps + alpha_eps * dealii::scalar_product(grad_u, grad_u))
                                  * fe_values.JxW(q);
          }

      }

  double bulk_energy = local_bulk_energy;
  double crack_energy = local_crack_energy;;

  std::cout << "Energies: " << timestep_number << "  " << time
        << "   " << bulk_energy
        << "   " << crack_energy
	    << std::endl;

}



