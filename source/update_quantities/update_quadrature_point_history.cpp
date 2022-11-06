#include "./../../include/dynamic_fracture.h"

  template <int dim>
  void Dynamic_Fracture_Problem<dim>::update_quadrature_point_history ()
  {
    QGauss<dim>   quadrature_formula(degree+2); 

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values | update_gradients | update_quadrature_points);
 
   std::vector<std::vector<Tensor<1,dim> > >
     old_solution_grads (quadrature_formula.size(),
			 std::vector<Tensor<1,dim> >(dim+1+dim));


   const Tensor<2, dim> Identity =
     ALE_Transformations::get_Identity<dim>();

   Tensor<2,dim> zero_matrix;
   zero_matrix.clear();



    for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
         cell != dof_handler.end(); ++cell)
      {
          PointHistory<dim> *local_quadrature_points_history
            = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
          Assert (local_quadrature_points_history >=
                  &  quadrature_point_history.front(),
                  ExcInternalError());
          Assert (local_quadrature_points_history <
                  &  quadrature_point_history.back(),
                  ExcInternalError());

          fe_values.reinit (cell);
          fe_values.get_function_gradients (solution,
                                            old_solution_grads);

          // Then loop over the quadrature points of this cell:
          for (unsigned int q=0; q<quadrature_formula.size(); ++q)
            {
              // On each quadrature point, compute the strain increment from
              // the gradients, and multiply it by the stress-strain tensor to
              // get the stress update. Then add this update to the already
              // existing strain at this point.

	     

	      const Tensor<2,dim> grad_u = ALE_Transformations 
		::get_grad_u<dim> (q, old_solution_grads);
	      

	      //std::cout << "Bin drin" << std::endl;
	      // Linearized strain
	      const Tensor<2,dim> E = 0.5 * (grad_u + transpose(grad_u));
	      
	      const double tr_E = Structure_Terms_in_ALE
		::get_tr_E<dim> (E);
	      
	      
	      Tensor<2,dim> stress_term;
	      stress_term.clear();
	      stress_term = lame_coefficient_lambda * tr_E * Identity
		+ 2 * lame_coefficient_mu * E;

	      Tensor<2,dim> stress_term_plus;
              Tensor<2,dim> stress_term_minus;
	      
	      if (timestep_number > 1)
		{
		  decompose_stress(stress_term_plus, stress_term_minus,
                                   E, tr_E, zero_matrix , 0.0,
                                   lame_coefficient_lambda,
                                   lame_coefficient_mu, false);
		}
	      else 
		{
		  stress_term_plus = stress_term;
		  stress_term_minus = 0;
		}
	      

	      if (bool_set_initial_strain_history && (test_case == "pressurized"))
		{
		  //std::cout << "Drin" << std::endl;
		  double top = 2.0 + alpha_eps/4.0;
		  double bottom = 2.0 - alpha_eps/4.0;
		  if (((fe_values.quadrature_point(q)[0] >= 1.8) && (fe_values.quadrature_point(q)[0] <= 2.2)) &&
		      ((fe_values.quadrature_point(q)[1] >= bottom) && (fe_values.quadrature_point(q)[1] <= top))
		      )
		    {
		      stress_term_plus[0][0] = 1.0e+10;
		      stress_term_plus[0][1] = 0.0;
		      stress_term_plus[1][0] = 0.0;
		      stress_term_plus[1][1] = 1.0e+10;
		    }
		  else 
		    {
		      stress_term_plus[0][0] = 0.0;
		      stress_term_plus[0][1] = 0.0;
		      stress_term_plus[1][0] = 0.0;
		      stress_term_plus[1][1] = 0.0;
		    }
		  
		}


	      if (bool_set_initial_strain_history && (test_case == "miehe_shear"))
		{
		  //std::cout << "Drin" << std::endl;
		   double top = 0.5 + alpha_eps/4.0;
		   double bottom = 0.5 - alpha_eps/4.0;
		 
		  if (((fe_values.quadrature_point(q)[0] >= 0.5) && (fe_values.quadrature_point(q)[0] <= 1.0)) &&
		      ((fe_values.quadrature_point(q)[1] >= bottom) && (fe_values.quadrature_point(q)[1]  <= top))
		      )
		    {
		      stress_term_plus[0][0] = 0.0; //1.0e+14;
		      stress_term_plus[0][1] = 0.0;
		      stress_term_plus[1][0] = 0.0;
		      stress_term_plus[1][1] = 0.0; //1.0e+14;
		    }
		  else 
		    {
		      stress_term_plus[0][0] = 0.0;
		      stress_term_plus[0][1] = 0.0;
		      stress_term_plus[1][0] = 0.0;
		      stress_term_plus[1][1] = 0.0;
		    }
		  
		}
	      	      


 
	      double stress_term_plus_norm = std::sqrt(stress_term_plus[0][0] * stress_term_plus[0][0] +
						       stress_term_plus[0][1] * stress_term_plus[0][1] +
						       stress_term_plus[1][0] * stress_term_plus[1][0] +
						       stress_term_plus[1][1] * stress_term_plus[1][1]);

	      Tensor<2,dim> old_stress = local_quadrature_points_history[q].old_stress;

	      double old_stress_norm = std::sqrt(old_stress[0][0] * old_stress[0][0] +
						 old_stress[0][1] * old_stress[0][1] +
						 old_stress[1][0] * old_stress[1][0] +
						 old_stress[1][1] * old_stress[1][1]);

	      Tensor<2,dim> new_stress;
	      if (stress_term_plus_norm > old_stress_norm)
		{
		  new_stress  = stress_term_plus; //old_stress + stress_term_plus;
		  //std::cout << "New stress" << std::endl;
		}
	      else
		{
		  new_stress  = old_stress;
		  //std::cout << "Oooooooooooooooold stress" << std::endl;
		}

              // The result of all these operations is then written back into
              // the original place:
              local_quadrature_points_history[q].old_stress
                = new_stress;
            }
        }
  }


