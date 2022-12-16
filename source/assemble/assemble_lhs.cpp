#include "./../../include/dynamic_fracture.h"

// In this function, we assemble the Jacobian matrix
// for the Newton iteration. 
template <int dim>
void Dynamic_Fracture_Problem<dim>::assemble_system_matrix ()
{
  #if DEAL_II_VERSION_GTE(9,0,0)
    timer.enter_subsection("Assemble Matrix.");
  #else
	timer.enter_section("Assemble Matrix.");
  #endif
  
  system_matrix=0;
     
  QGauss<dim>   quadrature_formula(degree+2);  
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values |
                           update_gradients);
  
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
				    update_values         | update_quadrature_points  |
				    update_normal_vectors | update_gradients |
				    update_JxW_values);
   
  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  
  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int n_face_q_points   = face_quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell); 
		

  // Now, we are going to use the 
  // FEValuesExtractors to determine
  // the four principle variables
  const FEValuesExtractors::Vector displacements (0); // 2
  const FEValuesExtractors::Scalar phase_field (dim); // 4
  const FEValuesExtractors::Vector velocities (dim+1); // 4

  // We declare Vectors and Tensors for 
  // the solutions at the previous Newton iteration:
  std::vector<Vector<double> > old_solution_values (n_q_points, 
				 		    Vector<double>(dim+1+dim));

  std::vector<std::vector<Tensor<1,dim> > > old_solution_grads (n_q_points, 
								std::vector<Tensor<1,dim> > (dim+1+dim));

  std::vector<Vector<double> >  old_solution_face_values (n_face_q_points, 
							  Vector<double>(dim+1+dim));
       
  std::vector<std::vector<Tensor<1,dim> > > old_solution_face_grads (n_face_q_points, 
								     std::vector<Tensor<1,dim> > (dim+1+dim));
    

  // We declare Vectors and Tensors for 
  // the solution at the previous time step:
   std::vector<Vector<double> > old_timestep_solution_values (n_q_points, 
				 		    Vector<double>(dim+1+dim));

   std::vector<Vector<double> >   old_old_timestep_solution_values (n_q_points, 
								    Vector<double>(dim+1+dim));


  std::vector<std::vector<Tensor<1,dim> > > old_timestep_solution_grads (n_q_points, 
  					  std::vector<Tensor<1,dim> > (dim+1+dim));


  std::vector<Vector<double> >   old_timestep_solution_face_values (n_face_q_points, 
								    Vector<double>(dim+1+dim));
  
    
  std::vector<std::vector<Tensor<1,dim> > >  old_timestep_solution_face_grads (n_face_q_points, 
									       std::vector<Tensor<1,dim> > (dim+1+dim));

std::vector<Vector<double> > old_solution_values_lambda_penal_func (n_q_points, 
								     Vector<double>(dim+1+dim));

   
  // Declaring test functions:
  std::vector<Tensor<1,dim> > phi_i_u (dofs_per_cell); 
  std::vector<Tensor<2,dim> > phi_i_grads_u(dofs_per_cell);
  std::vector<double>         phi_i_pf(dofs_per_cell); 
  std::vector<Tensor<1,dim> > phi_i_grads_pf (dofs_per_cell);
  // u and v here components of displacements
  std::vector<Tensor<1,dim> > phi_i_v (dofs_per_cell); 
  std::vector<Tensor<2,dim> > phi_i_grads_v(dofs_per_cell);

  // This is the identity matrix in two dimensions:
  const Tensor<2,dim> Identity = ALE_Transformations
    ::get_Identity<dim> ();

  Tensor<2,dim> zero_matrix;
  zero_matrix.clear();
 				     				   
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  unsigned int cell_counter = 0;
  for (; cell!=endc; ++cell)
    { 
      //lame_coefficient_mu = lame_coefficient_mu_vector(cell_counter);
      //lame_coefficient_lambda = lame_coefficient_lambda_vector(cell_counter);
      cell_counter ++;

      fe_values.reinit (cell);
      local_matrix = 0;
      
      // We need the cell diameter to control the fluid mesh motion
      cell_diameter = cell->diameter();

      if (current_test_case.Is_het3d())
	{
	  E_modulus = func_emodulus->value(cell->center(), 0);
	  E_modulus += 1.0;
	  //E_modulus = 5.0;
	  
	  lame_coefficient_mu = E_modulus / (2.0 * (1 + poisson_ratio_nu));
	  
	  lame_coefficient_lambda = (2 * poisson_ratio_nu * lame_coefficient_mu)
	    / (1.0 - 2 * poisson_ratio_nu);
	}

      // Old Newton iteration values
      fe_values.get_function_values (solution, old_solution_values);
      fe_values.get_function_gradients (solution, old_solution_grads);
      
      // Old_timestep_solution values
      fe_values.get_function_values (old_timestep_solution, old_timestep_solution_values);
      fe_values.get_function_gradients (old_timestep_solution, old_timestep_solution_grads);

      // Old_Old_timestep_solution values
      fe_values.get_function_values (old_old_timestep_solution, old_old_timestep_solution_values);
  
   // Old Newton iteration values lambda 
      fe_values.get_function_values (solution_lambda_penal_func, old_solution_values_lambda_penal_func);
   
  
  const PointHistory<dim> *local_quadrature_points_data
	= reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

	  for (unsigned int q=0; q<n_q_points; ++q)
	    {	      
	      for (unsigned int k=0; k<dofs_per_cell; ++k)
		{
		  phi_i_u[k]       = fe_values[displacements].value (k, q);
		  phi_i_grads_u[k] = fe_values[displacements].gradient (k, q);
		  phi_i_pf[k]       = fe_values[phase_field].value (k, q);
		  phi_i_grads_pf[k] = fe_values[phase_field].gradient (k, q);
		  phi_i_v[k]       = fe_values[velocities].value (k, q);
		  phi_i_grads_v[k] = fe_values[velocities].gradient (k, q);
		}
	      

	      const double pf = old_solution_values[q](dim); 
	      const double old_timestep_pf = old_timestep_solution_values[q](dim); 
	      const double old_old_timestep_pf = old_old_timestep_solution_values[q](dim); 

	      const double  lambda_penal_func = old_solution_values_lambda_penal_func[q](dim);

	      double pf_extra = pf;
              // Linearization by extrapolation to cope with non-convexity of the underlying
              // energy functional.
              // This idea might be refined in a future work (be also careful because
              // theoretically, we do not have time regularity; therefore extrapolation in time
              // might be questionable. But for the time being, this is numerically robust.
              pf_extra = old_old_timestep_pf + (time - (time-old_timestep-old_old_timestep))/
                         (time-old_timestep - (time-old_timestep-old_old_timestep)) * (old_timestep_pf - old_old_timestep_pf);


              if (pf_extra <= 0.0)
                pf_extra = 0.0;
              if (pf_extra >= 1.0)
                pf_extra = 1.0;

	      const Tensor<2,dim> grad_u = ALE_Transformations 
		::get_grad_u<dim> (q, old_solution_grads);

	      double divergence_u = old_solution_grads[q][0][0] +  old_solution_grads[q][1][1];
	      if (dim == 3)
		divergence_u += old_solution_grads[q][2][2];

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
	      
	      if ((timestep_number > 1) && (dim == 2) && bool_use_stress_splitting)
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

	      const Tensor<2,dim> &stress_term_history
		= local_quadrature_points_data[q].old_stress;

	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  // Simple penalization
//		  double pf_minus_old_timestep_pf_plus = 0.0;
//		  if ((pf - old_timestep_pf) < 0.0)
//		    pf_minus_old_timestep_pf_plus = 0.0;
//		  else 
//		    pf_minus_old_timestep_pf_plus = phi_i_pf[i]; 


		  double chi = 0.0;
		  if ((lambda_penal_func + gamma_penal * (pf - old_timestep_pf)) > 0.0)
		    chi = 1.0;
		  else 
		    chi = 0.0;

		  double divergence_u_LinU = phi_i_grads_u[i][0][0] + phi_i_grads_u[i][1][1];
		  if (dim == 3)
		    divergence_u_LinU += phi_i_grads_u[i][2][2];
	    	     
		
		  const Tensor<2, dim> E_LinU = 0.5
		    * (phi_i_grads_u[i] + transpose(phi_i_grads_u[i]));

		  const double tr_E_LinU = dealii::trace(E_LinU);

		  Tensor<2,dim> stress_term_LinU;
		  stress_term_LinU = lame_coefficient_lambda * tr_E_LinU * Identity
		    + 2 * lame_coefficient_mu * E_LinU;

		  Tensor<2,dim> stress_term_plus_LinU;
		  Tensor<2,dim> stress_term_minus_LinU;
		  if ((timestep_number > 1) && (dim == 2) && bool_use_stress_splitting)
		    {
		  decompose_stress(stress_term_plus_LinU, stress_term_minus_LinU,
				   E, tr_E, E_LinU, tr_E_LinU,
				   lame_coefficient_lambda,
				   lame_coefficient_mu,
				   true);
		    }
		  else 
		    {
		      stress_term_plus_LinU = stress_term_LinU;
		      stress_term_minus_LinU = 0;
		    }


		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    {
		      // STVK 
		      const unsigned int comp_j = fe.system_to_component_index(j).first; 
		      if (comp_j < dim)
			{
			  if (bool_use_dynamic_code_with_velocities)
			    {
			      // TODO: density, dividing over time step
			      local_matrix(j,i) += density_structure * phi_i_v[i] * phi_i_u[j] * fe_values.JxW(q);
			    }

			  // pf is solution variable
			  if (!bool_use_pf_extra)
			    {
			      // Simple Hessian modification (not very reliable)
			      //local_matrix(j,i) += 0.0 * timestep * theta * 
			      //	1.0/(cell_diameter*cell_diameter) * 0.1 * lame_coefficient_mu * (phi_i_u[i] * phi_i_u[j]
			      //	 ) * fe_values.JxW(q);



			      local_matrix(j,i) += timestep * theta * 
				(delta_fixed_point_newton * dealii::scalar_product((1-constant_k) * 2.0 * pf * phi_i_pf[i] * stress_term_plus,phi_i_grads_u[j])
				 + dealii::scalar_product(((1-constant_k) * pf * pf + constant_k) * stress_term_plus_LinU, phi_i_grads_u[j])
				 + dealii::scalar_product(stress_term_minus_LinU, phi_i_grads_u[j])
				 // Pressure (pf is solution variable)
				 - 2.0 * (alpha_biot - 1.0) * current_pressure * pf * phi_i_pf[i] * divergence_u_LinU
				 ) * fe_values.JxW(q);
			    }
			  else if (bool_use_pf_extra)
			    {
			      // pf extrapolated
			      local_matrix(j,i) += timestep * theta * 
				(dealii::scalar_product(((1-constant_k) * pf_extra * pf_extra + constant_k) * 	      
						stress_term_plus_LinU, phi_i_grads_u[j])
				 + dealii::scalar_product(stress_term_minus_LinU, phi_i_grads_u[j])
				 ) * fe_values.JxW(q);
			    }
			}		     
		      else if (comp_j == dim)
			{

			  if (!bool_use_strain_history)
			    {

			  // Simple penalization
			  //local_matrix(j,i) += delta_penal *  1.0/cell_diameter *  phi_i_pf[j] 
			  //  * pf_minus_old_timestep_pf_plus* fe_values.JxW(q);

			  // Augmented Lagrangian penalization
			      local_matrix(j,i) += chi * gamma_penal * phi_i_pf[i] * phi_i_pf[j] * fe_values.JxW(q);


			      local_matrix(j,i) += timestep * theta *
				( (1-constant_k) * (dealii::scalar_product(stress_term_plus_LinU, E)
						    + dealii::scalar_product(stress_term_plus, E_LinU)) * pf * phi_i_pf[j]
				  +(1-constant_k) * dealii::scalar_product(stress_term_plus, E) * phi_i_pf[i] * phi_i_pf[j]
				  + G_c/alpha_eps * phi_i_pf[i] * phi_i_pf[j]  
				  + G_c * alpha_eps * phi_i_grads_pf[i] * phi_i_grads_pf[j]
				  // Pressure terms
				  - 2.0 * (alpha_biot - 1.0) * current_pressure *
				  (pf * divergence_u_LinU + phi_i_pf[i] * divergence_u) * phi_i_pf[j]
				  ) * fe_values.JxW(q);      

			    }
			  else if (bool_use_strain_history)
			    {

			      local_matrix(j,i) += timestep * theta *
			    ((1-constant_k) * (0.0 * dealii::scalar_product(stress_term_plus_LinU, E)
			     		       + dealii::scalar_product(stress_term_history, E_LinU)) * pf * phi_i_pf[j]
			     +(1-constant_k) * dealii::scalar_product(stress_term_history, E) * phi_i_pf[i] * phi_i_pf[j]

			     + G_c/alpha_eps * phi_i_pf[i] * phi_i_pf[j]  
			     + G_c * alpha_eps * phi_i_grads_pf[i] * phi_i_grads_pf[j]
			     // Pressure terms
			     - 2.0 * (alpha_biot - 1.0) * current_pressure *
			     (pf * divergence_u_LinU + phi_i_pf[i] * divergence_u) * phi_i_pf[j]
			     ) * fe_values.JxW(q);     

			    }

			}
		      else if (comp_j > dim)
			{
			  //local_matrix(j,i) += scalar_product(phi_i_grads_v[i],phi_i_grads_v[j])
			  // * fe_values.JxW(q);

			  if (bool_use_dynamic_code_with_velocities)
			    {
			      local_matrix(j,i) += (phi_i_u[i] - timestep * theta * phi_i_v[i])  * phi_i_v[j] * fe_values.JxW(q);
			    }
			  else 
			    {
			      local_matrix(j,i) += phi_i_v[i]  * phi_i_v[j]
				* fe_values.JxW(q);
			    }


			}

		      
		    }  // end j dofs
		  	     
		}  // end i dofs 
	      
	    }   // end n_q_points  

	  
	  cell->get_dof_indices (local_dof_indices);
	  constraints.distribute_local_to_global (local_matrix, local_dof_indices,
						  system_matrix);

      // end cell
    }   
  
  #if DEAL_II_VERSION_GTE(9,0,0)
    timer.leave_subsection("Assemble Matrix.");
  #else
	timer.exit_section();
  #endif
  
}