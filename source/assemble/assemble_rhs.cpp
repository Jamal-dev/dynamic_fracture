#include "./../../include/dynamic_fracture.h"

// In this function we assemble the semi-linear 
// of the right hand side of Newton's method (its residual).
// The framework is in principal the same as for the 
// system matrix.
template <int dim>
void
Dynamic_Fracture_Problem<dim>::assemble_system_rhs ()
{
  timer.enter_section("Assemble Rhs.");
  system_rhs=0;
  
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
 
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
  const FEValuesExtractors::Vector displacements (0); 
  const FEValuesExtractors::Scalar phase_field (dim);
  const FEValuesExtractors::Vector velocities (dim+1);
 
  // old Newton step values
  std::vector<Tensor<1, dim> > old_displacement_values(n_q_points);
  std::vector<Tensor<1, dim> > old_velocity_values(n_q_points);
  std::vector<double> old_phase_field_values(n_q_points);

  // old Newton step grads
  std::vector<Tensor<2, dim> > old_displacement_grads(n_q_points);
  std::vector<Tensor<2, dim> > old_velocity_grads(n_q_points);
  std::vector<Tensor<1,dim> > old_phase_field_grads(n_q_points);

  /*old time step */
  std::vector<Tensor<1, dim> > old_timestep_displacement_values(n_q_points);
  std::vector<Tensor<1, dim> > old_timestep_velocity_values(n_q_points);
  std::vector<double> old_timestep_phase_field_values(n_q_points);

  // Old time step grads
  std::vector<Tensor<2,dim> > old_timestep_displacement_grads (n_q_points);
  std::vector<Tensor<2,dim> > old_timestep_velocity_grads (n_q_points);
  std::vector<Tensor<1,dim> > old_timestep_phase_field_grads (n_q_points);

  // old old time step values
  std::vector<Tensor<1, dim> > old_old_timestep_displacement_values(n_q_points);
  std::vector<Tensor<1, dim> > old_old_timestep_velocity_values(n_q_points);
  std::vector<double> old_old_timestep_phase_field_values(n_q_points);
  
  std::vector<Vector<double> > 
    old_solution_values (n_q_points, Vector<double>(dim+1+dim));

  std::vector<std::vector<Tensor<1,dim> > > 
    old_solution_grads (n_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));


  std::vector<Vector<double> > 
    old_solution_face_values (n_face_q_points, Vector<double>(dim+1+dim));
  
  std::vector<std::vector<Tensor<1,dim> > > 
    old_solution_face_grads (n_face_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));
  
  std::vector<Vector<double> > 
    old_timestep_solution_values (n_q_points, Vector<double>(dim+1+dim));

  std::vector<Vector<double> > 
    old_old_timestep_solution_values (n_q_points, Vector<double>(dim+1+dim));

  std::vector<std::vector<Tensor<1,dim> > > 
    old_timestep_solution_grads (n_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));

  std::vector<Vector<double> > 
    old_timestep_solution_face_values (n_face_q_points, Vector<double>(dim+1+dim));
     
  std::vector<std::vector<Tensor<1,dim> > > 
    old_timestep_solution_face_grads (n_face_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));

 std::vector<Vector<double> > old_solution_values_lambda_penal_func (n_q_points, 
								     Vector<double>(dim+1+dim));
 
   
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
      //Required for heterogeneous media
      //lame_coefficient_mu = lame_coefficient_mu_vector(cell_counter);
      //lame_coefficient_lambda = lame_coefficient_lambda_vector(cell_counter);
      cell_counter++;


      fe_values.reinit (cell);	 
      local_rhs = 0;   	
      
      cell_diameter = cell->diameter();

      if (test_case == "Het3D")
	{
	  E_modulus = func_emodulus->value(cell->center(), 0);
	  E_modulus += 1.0;
	  //E_modulus = 5.0;
	  
	  lame_coefficient_mu = E_modulus / (2.0 * (1 + poisson_ratio_nu));
	  
	  lame_coefficient_lambda = (2 * poisson_ratio_nu * lame_coefficient_mu)
	    / (1.0 - 2 * poisson_ratio_nu);
	}

      // old Newton iteration values
	  fe_values[displacements].get_function_values (solution, old_displacement_values);
	  fe_values[phase_field].get_function_values (solution, old_phase_field_values);
	  fe_values[velocities].get_function_values (solution, old_velocity_values);
	  
	  // old Newton iteration gradients
	  fe_values[displacements].get_function_gradients (solution, old_displacement_grads);
	  fe_values[phase_field].get_function_gradients (solution, old_phase_field_grads);
	  fe_values[velocities].get_function_gradients (solution, old_velocity_grads);

	  // old time step values
	  fe_values[displacements].get_function_values (old_timestep_solution, old_timestep_displacement_values);
	  fe_values[phase_field].get_function_values (old_timestep_solution, old_timestep_phase_field_values);
	  fe_values[velocities].get_function_values (old_timestep_solution, old_timestep_velocity_values);

	  // old time step gradients
	  fe_values[displacements].get_function_gradients (old_timestep_solution, old_timestep_displacement_grads);
	  fe_values[phase_field].get_function_gradients (old_timestep_solution, old_timestep_phase_field_grads);
	  fe_values[velocities].get_function_gradients (old_timestep_solution, old_timestep_velocity_grads);

	  // old old time step values
	  fe_values[displacements].get_function_values (old_old_timestep_solution, old_old_timestep_displacement_values);
	  fe_values[phase_field].get_function_values (old_old_timestep_solution, old_old_timestep_phase_field_values);
	  fe_values[velocities].get_function_values (old_old_timestep_solution, old_old_timestep_velocity_values);
      // old Newton iteration
      fe_values.get_function_values (solution, old_solution_values);
      fe_values.get_function_gradients (solution, old_solution_grads);
            
      // old timestep iteration
      fe_values.get_function_values (old_timestep_solution, old_timestep_solution_values);
      fe_values.get_function_gradients (old_timestep_solution, old_timestep_solution_grads);

      // old timestep iteration
      fe_values.get_function_values (old_old_timestep_solution, old_old_timestep_solution_values);

      // Old Newton iteration values lambda 
      fe_values.get_function_values (solution_lambda_penal_func, old_solution_values_lambda_penal_func);



    const PointHistory<dim> *local_quadrature_points_data
			= reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

      

	  for (unsigned int q=0; q<n_q_points; ++q)
	    {
	      const double pf = old_phase_field_values[q];
	      const double old_timestep_pf = old_timestep_phase_field_values[q];
		
		  const double old_old_timestep_pf = old_old_timestep_phase_field_values[q];
	      
	      
	      const double  lambda_penal_func = old_solution_values_lambda_penal_func[q](dim);

	      double pf_extra = pf;
              // Linearization by extrapolation to cope with non-convexity of the underlying
              // energy functional.
              // This idea might be refined in a future work (be also careful because
              // theoretically, we do not have time regularity; therefore extrapolation in time
              // might be questionable. But for the time being, this is numerically robust.
              pf_extra = old_old_timestep_pf + 
			  			(time - (time-old_timestep-old_old_timestep))/
                         (time-old_timestep - (time-old_timestep-old_old_timestep)) 
						 * (old_timestep_pf - old_old_timestep_pf);


              if (pf_extra <= 0.0)
                pf_extra = 0.0;
              if (pf_extra >= 1.0)
                pf_extra = 1.0;


	      // Simple penalization
//	      double pf_minus_old_timestep_pf_plus = 0.0;
//	      if ((pf - old_timestep_pf) < 0.0)
//		pf_minus_old_timestep_pf_plus = 0.0;
//	      else 
//		pf_minus_old_timestep_pf_plus = pf - old_timestep_pf;


	      // Augmented Lagrangian; but first term missing
	      double pf_minus_old_timestep_pf_plus = 0.0;
	      if ((lambda_penal_func + gamma_penal * (pf - old_timestep_pf)) < 0.0)
		pf_minus_old_timestep_pf_plus = 0.0;
	      else 
		pf_minus_old_timestep_pf_plus = lambda_penal_func + gamma_penal * (pf - old_timestep_pf);

	      //std::cout << pf_minus_old_timestep_pf_plus << std::endl;

	    //   const Tensor<1,dim> grad_pf = ALE_Transformations
		// ::get_grad_pf<dim> (q, old_solution_grads);
		const Tensor<1,dim> grad_pf = old_phase_field_grads[q];

	    //   const Tensor<1,dim> old_timestep_grad_pf = ALE_Transformations
		// ::get_grad_pf<dim> (q, old_timestep_solution_grads);

		const Tensor<1,dim> old_timestep_grad_pf = 
									old_timestep_phase_field_grads[q];

	    //   const Tensor<1,dim> v = ALE_Transformations 
		// ::get_v<dim> (q, old_solution_values);
		const Tensor<1,dim> v = old_velocity_values[q];

	    //   const Tensor<1,dim> old_timestep_v = ALE_Transformations 
		// ::get_v<dim> (q, old_timestep_solution_values);

		const Tensor<1,dim> old_timestep_v = old_timestep_velocity_values[q];

	    //   const Tensor<2,dim> grad_v = ALE_Transformations 
		// ::get_grad_v<dim> (q, old_solution_grads);

		const Tensor<2,dim> grad_v = old_velocity_grads[q];
	      

	    //   const Tensor<1,dim> u = ALE_Transformations 
		// ::get_u<dim> (q, old_solution_values);

		const Tensor<1,dim> u = old_displacement_values[q];
	      
	    //   const Tensor<1,dim> old_timestep_u = ALE_Transformations 
		// ::get_u<dim> (q, old_timestep_solution_values);

		const Tensor<1,dim> old_timestep_u = old_timestep_displacement_values[q];
	      
	    //   const Tensor<2,dim> grad_u = ALE_Transformations 
		// ::get_grad_u<dim> (q, old_solution_grads);

		const Tensor<2,dim> grad_u = old_displacement_grads[q];

	    //   const Tensor<2,dim> old_timestep_grad_u = ALE_Transformations 
		// ::get_grad_u<dim> (q, old_timestep_solution_grads);

		const Tensor<2,dim> old_timestep_grad_u = old_timestep_displacement_grads[q];

		double divergence_u = ALE_Transformations::get_divergence_u<dim> (grad_u);	      
	    //   double divergence_u = old_solution_grads[q][0][0] +  old_solution_grads[q][1][1];
	    //   if (dim == 3)
		// divergence_u += old_solution_grads[q][2][2];

	    double old_timestep_divergence_u = ALE_Transformations::get_divergence_u<dim> (old_timestep_grad_u);
		//   double old_timestep_divergence_u = old_timestep_solution_grads[q][0][0] +  old_timestep_solution_grads[q][1][1];
	    //   if (dim == 3)
		// old_timestep_divergence_u += old_timestep_solution_grads[q][2][2];

	      // Linearized strain
	      const Tensor<2,dim> E = 0.5 * (grad_u + transpose(grad_u));
	      
	      const double tr_E = Structure_Terms_in_ALE
								::get_tr_E<dim> (E);

	      const Tensor<2,dim> old_timestep_E = 0.5 * (old_timestep_grad_u + transpose(old_timestep_grad_u));
	      
	      const double old_timestep_tr_E = Structure_Terms_in_ALE
										::get_tr_E<dim> (old_timestep_E);
	      
	      
	      Tensor<2,dim> stress_term;
	      stress_term.clear();
	      stress_term = lame_coefficient_lambda * tr_E * Identity
							+ 2 * lame_coefficient_mu * E;

	      Tensor<2,dim> old_timestep_stress_term;
	      old_timestep_stress_term.clear();
	      old_timestep_stress_term = lame_coefficient_lambda * old_timestep_tr_E * Identity
								+ 2 * lame_coefficient_mu * old_timestep_E;

	      Tensor<2,dim> stress_term_plus;
              Tensor<2,dim> stress_term_minus;

	      Tensor<2,dim> old_timestep_stress_term_plus;
              Tensor<2,dim> old_timestep_stress_term_minus;
	      
	      if ((timestep_number > 1) && (dim == 2) && bool_use_stress_splitting)
		{
		  decompose_stress(stress_term_plus, stress_term_minus,
                                   E, tr_E, zero_matrix , 0.0,
                                   lame_coefficient_lambda,
                                   lame_coefficient_mu, false);

		  // Old timestep
		  decompose_stress(old_timestep_stress_term_plus, old_timestep_stress_term_minus,
                                   old_timestep_E, old_timestep_tr_E, zero_matrix , 0.0,
                                   lame_coefficient_lambda,
                                   lame_coefficient_mu, false);
		}
	      else 
		{
		  stress_term_plus = stress_term;
		  stress_term_minus = 0;

		  // Old timestep
		  old_timestep_stress_term_plus = old_timestep_stress_term;
		  old_timestep_stress_term_minus = 0;
		}
	      

	       const Tensor<2,dim> &stress_term_history
					= local_quadrature_points_data[q].old_stress; 


	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  // Info: Compute velocities
		  const unsigned int comp_i = fe.system_to_component_index(i).first; 
		  if (comp_i < dim) // displacement part
		    { 
		      const Tensor<1,dim> phi_i_u = fe_values[displacements].value (i, q);
		      const Tensor<2,dim> phi_i_grads_u = fe_values[displacements].gradient (i, q);

		      double divergence_u_LinU = phi_i_grads_u[0][0] + phi_i_grads_u[1][1];
		      if (dim == 3)
			divergence_u_LinU += phi_i_grads_u[2][2];


		      // Acceleration term
		      if (bool_use_dynamic_code_with_velocities)
			{
			//   local_rhs(i) -= (v - old_timestep_v) * phi_i_u * fe_values.JxW(q);
			  local_rhs(i) += density_structure * (old_timestep_v - v) * phi_i_u * fe_values.JxW(q);
			}

		      // pf is solution variable
		      if (!bool_use_pf_extra)
			{
			  
			  // Current timestep solution
			  local_rhs(i) -= timestep * theta * 
			    (dealii::scalar_product(((1-constant_k) * pf * pf + constant_k) *	  
					    stress_term_plus, phi_i_grads_u)
			     + dealii::scalar_product(stress_term_minus, phi_i_grads_u)
			     // Pressure terms (pf is solution variable)
			     - (alpha_biot - 1.0) * current_pressure * pf * pf * divergence_u_LinU
			     ) * fe_values.JxW(q); 

			  // Old timestep solution
			  local_rhs(i) -= timestep * (1.0-theta) * 
			    (dealii::scalar_product(((1-constant_k) * old_timestep_pf * old_timestep_pf + constant_k) *	  
					    old_timestep_stress_term_plus, phi_i_grads_u)
			     + dealii::scalar_product(old_timestep_stress_term_minus, phi_i_grads_u)
			     // Pressure terms (pf is solution variable)
			     - (alpha_biot - 1.0) * old_timestep_current_pressure * old_timestep_pf * old_timestep_pf * divergence_u_LinU
			     ) * fe_values.JxW(q); 

			  



			}
		      else if (bool_use_pf_extra)
			{
			 
			  // Current timestep solution
			  // pf extrapolated
			  local_rhs(i) -= timestep * theta * 
			    (dealii::scalar_product(((1-constant_k) * pf_extra * pf_extra + constant_k) *	  
					    stress_term_plus, phi_i_grads_u)
			     + dealii::scalar_product(stress_term_minus, phi_i_grads_u)
			     // Pressure terms (extrapolated)
			     - (alpha_biot - 1.0) * current_pressure * pf_extra * pf_extra * divergence_u_LinU
			     ) * fe_values.JxW(q); 

			  // Old timestep solution
			  // TODO: need to define previous timestep extrapolation
			  local_rhs(i) -= timestep * (1.0 - theta) * 
			    (dealii::scalar_product(
							(   
								(1-constant_k) * pf_extra * pf_extra 
								+ constant_k) 
								* old_timestep_stress_term_plus, 
											phi_i_grads_u)
			     + dealii::scalar_product(old_timestep_stress_term_minus
				 							, phi_i_grads_u)
			     // Pressure terms (extrapolated)
			     - (alpha_biot - 1.0) * old_timestep_current_pressure * pf_extra * pf_extra * divergence_u_LinU
			     ) * fe_values.JxW(q); 
			  
			  
			} 
		      
		    }		
		  else if (comp_i == dim)
		    {
		      const double phi_i_pf = fe_values[phase_field].value (i, q);
		      const Tensor<1,dim> phi_i_grads_pf = fe_values[phase_field].gradient (i, q);

		      if (!bool_use_strain_history)
			{
			  //  Simple penalization
			//   local_rhs(i) -= delta_penal *  1.0/(cell_diameter * cell_diameter)  
			//   	* pf_minus_old_timestep_pf_plus * phi_i_pf * fe_values.JxW(q);

		     
			  // Augmented Lagrangian penalization
			  local_rhs(i) -= pf_minus_old_timestep_pf_plus * phi_i_pf * fe_values.JxW(q);
		      
			  // Current time step
			  local_rhs(i) -= timestep * theta * 
			    ((1.0 - constant_k) * dealii::scalar_product(stress_term_plus, E) * pf * phi_i_pf
			     - G_c/alpha_eps * (1.0 - pf) * phi_i_pf
			     + G_c * alpha_eps * grad_pf * phi_i_grads_pf
			     // Pressure terms
			     - 2.0 * (alpha_biot - 1.0) * current_pressure * pf * divergence_u * phi_i_pf
			     ) * fe_values.JxW(q);

			  // Old timestep
			  local_rhs(i) -= timestep * (1.0-theta) * 
			    ((1.0 - constant_k) 
				* dealii::scalar_product(old_timestep_stress_term_plus, E) 
				* old_timestep_pf * phi_i_pf
			     - G_c/alpha_eps * (1.0 - old_timestep_pf) * phi_i_pf
			     + G_c * alpha_eps * old_timestep_grad_pf * phi_i_grads_pf
			     // Pressure terms
			     - 2.0 * (alpha_biot - 1.0) * old_timestep_current_pressure * old_timestep_pf * old_timestep_divergence_u * phi_i_pf
			     ) * fe_values.JxW(q);
			  
			}
		      else if (bool_use_strain_history)
			{
			  std::cout << "Aborting ..." << std::endl;
			  abort();
			  // only implemented for BE (backward Euler) so far
//                         local_rhs(i) -= timestep * theta * 
//			((1.0 - constant_k) * scalar_product(stress_term_history, E) * pf * phi_i_pf
//			 - G_c/alpha_eps * (1.0 - pf) * phi_i_pf
//			 + G_c * alpha_eps * grad_pf * phi_i_grads_pf
//			 // Pressure terms
//			 - 2.0 * (alpha_biot - 1.0) * current_pressure * pf * divergence_u * phi_i_pf
//			 ) * fe_values.JxW(q);

			}

		    } 
		  else if (comp_i > dim)
		    { 
		      // Info: Compute displacements
		      const Tensor<1,dim> phi_i_v       = fe_values[velocities].value (i, q);
		      const Tensor<2,dim> phi_i_grads_v = fe_values[velocities].gradient (i, q);
		      
		      //local_rhs(i) -= 
		      //	(scalar_product(grad_v,phi_i_grads_v) - 1.0 * phi_i_v[0]
		      //	 ) * fe_values.JxW(q); 

		      if (bool_use_dynamic_code_with_velocities)
			{
			  //local_rhs(i) -= ((u - old_timestep_u) - timestep * (theta * v + (1-theta) * old_timestep_v)) * phi_i_v * fe_values.JxW(q);
			  // desnity is common so not using it; it should not be used in
			  // left side as well
			  local_rhs(i) += (
									(old_timestep_u - u)  
			  					+ timestep *(theta * v + (1-theta) * old_timestep_v))
								 * phi_i_v * fe_values.JxW(q);  
			}
		      else 
			{
			  local_rhs(i) -= 
			    (v * phi_i_v
			     ) * fe_values.JxW(q); 
			}

		    }  

		  	  
		} // end i
	      		   
	    } // end n_q_points 
	  


	  // upper traction (non-pay zone)
	  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	    {
	      if (cell->face(face)->at_boundary() && 		  
		  (cell->face(face)->boundary_id() == 3) 
		  )
		{
		  
		  fe_face_values.reinit (cell, face);
		  
		  for (unsigned int q=0; q<n_face_q_points; ++q)
		    {	
		      Tensor<1,dim> neumann_value;
		      neumann_value[0] = 0.0; //time * traction_x;
		      neumann_value[1] = 0.0; //time * traction_y;
			
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
			  const unsigned int comp_i = fe.system_to_component_index(i).first; 
			  if (comp_i == 0 || comp_i == 1)
			    {  
			      local_rhs(i) +=  1.0 * (timestep * theta * 
						      neumann_value * fe_face_values[displacements].value (i, q) 
						      ) * fe_face_values.JxW(q);					   
			    }
			  // end i
			}  
		      // end face_n_q_points    
		    }                                     
		} 
	    }  // end face integrals upper fraction



	  // lower traction (non-pay zone)
	  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	    {
	      if (cell->face(face)->at_boundary() && 		  
		  (cell->face(face)->boundary_id() == 2) 
		  )
		{
		  
		  fe_face_values.reinit (cell, face);
		  
		  for (unsigned int q=0; q<n_face_q_points; ++q)
		    {	
		      Tensor<1,dim> neumann_value;
		      neumann_value[0] = - 0.0; //time * traction_x;
		      neumann_value[1] = - 0.0; //time * traction_y;
			
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
			  const unsigned int comp_i = fe.system_to_component_index(i).first; 
			  if (comp_i == 0 || comp_i == 1)
			    {  
			      local_rhs(i) +=  1.0 * (timestep * theta * 
						      neumann_value * fe_face_values[displacements].value (i, q) 
						      ) * fe_face_values.JxW(q);					   
			    }
			  // end i
			}  
		      // end face_n_q_points    
		    }                                     
		} 
	    }  // end face integrals lower traction










	  cell->get_dof_indices (local_dof_indices);
	  constraints.distribute_local_to_global (local_rhs, local_dof_indices,
						  system_rhs);
	  

      
    }  // end cell
      
  timer.exit_section();
}


