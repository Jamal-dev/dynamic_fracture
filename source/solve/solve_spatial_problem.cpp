#include "./../../include/dynamic_fracture.h"

template <int dim>
void Dynamic_Fracture_Problem<dim>::solve_spatial_problem () 
{ 

  BlockVector<double> tmp_old_sol;
  tmp_old_sol = solution;

  unsigned int prec_corr_counter_output_a = 100000 + timestep_number;
  unsigned int prec_corr_counter_output_b = 200000 + timestep_number;
  unsigned int prec_corr_counter_output_c = 300000 + timestep_number;

    redo_step:
      std::cout << "Timestep " << timestep_number 
		<< " (" << time_stepping_scheme 
		<< ")" <<    ": " << time
		<< " (" << timestep << ")"
		<< "   " << "Cells: " << triangulation.n_active_cells()
		<< "   " << "DoFs: " << dof_handler.n_dofs()
		<< "\n======================================" 
		<< "==========================================" 
		<< std::endl;
      if (!bool_use_error_oriented_Newton)
	{
	  std::cout << "Iter\t" << "Residual\t"
		    << "Red. New/Old\t" 
		    << "Bui.Mat\t"
		    << "LS\t"
		    << "Delta (Fix-P)\t"
		    << "CPU time"
		    << std::endl;
	}
      else if (bool_use_error_oriented_Newton)
	{
	  std::cout << "Iter\t" << "Update\t\t"
		    << "Update(bar)\t" 
		    << "Residual\t"
		    << "Red. Up/Old-Up\t"
		    << "Bui.Mat\t"
		    << "CPU time"
		    << std::endl;


	}


      BlockVector<double> solution_lambda_penal_func_difference;
      penal_iterations = 0;
      double L2_error = 0.0;
      double L2_error_relative = 0.0;
      double initial_L2_norm = 0.0;

      /*
	// TODO
      if (test_case == "miehe_shear" ||
	  test_case == "pressurized" ||
	  test_case == "screw_domi" ||
	  test_case == "Sneddon")
	{
	  lower_bound_newton_update = 1.0e-8;
	  lower_bound_newton_residuum = 1.0e-8; 
	}
      else if (test_case == "l_shaped")
	{
	  lower_bound_newton_update = 1.0e-8;
	  lower_bound_newton_residuum = 1.0e-8; 
	}
      */

      double pressure_increment_factor = 0.1;
      
 pf_extra_step:
      
      // Reset penal function for new time step
      solution_lambda_penal_func = 0;

      // Augmented Lagrangian iteration
      do
	{
	  // Step 1: Solve forward system
	   set_initial_guess_Newton (time);
	  if (!bool_use_error_oriented_Newton)
	    {
	      delta_fixed_point_newton = 1.0;
	      newton_iteration (time);
	    }
	  else if (bool_use_error_oriented_Newton)
	    {
	      // Only for initialization - apart from this,
	      // delta_fixed_point Newton is not used
	      // in the error oriented Newton method
	      delta_fixed_point_newton = 1.0;
	      newton_iteration_error_based (time);   
	    }
 

	  
	  // Step 2: Update lambda
	  old_timestep_solution_lambda_penal_func = solution_lambda_penal_func;
	  for (unsigned int j=0; j<solution_lambda_penal_func.size(); ++j)
	    {
	      if ((solution_lambda_penal_func(j)  + gamma_penal * (solution(j)  - old_timestep_solution(j) )) > 0.0)
		solution_lambda_penal_func(j) 
		  = solution_lambda_penal_func(j)  + gamma_penal * (solution(j)  - old_timestep_solution(j) );
	      else 
		solution_lambda_penal_func(j) = 0.0;

	      //std::cout << solution_lambda_penal_func(j) << std::endl;
	    }

	  // Step 3: Check stopping criterium
	  solution_lambda_penal_func_difference = solution_lambda_penal_func;
	  solution_lambda_penal_func_difference -= old_timestep_solution_lambda_penal_func;
	  

	  ComponentSelectFunction<dim> value_select (dim, dim+1+dim);
	  Vector<float> difference_per_cell (triangulation.n_active_cells());
	  VectorTools::integrate_difference (dof_handler,
					     //tmp_vec,
					     solution_lambda_penal_func_difference,
					     ZeroFunction<dim>(dim+1+dim),					     
					     difference_per_cell,
					     QGauss<dim>(fe.degree+1),
					     VectorTools::L2_norm,
					     &value_select);
	  L2_error = difference_per_cell.l2_norm();

	  if (penal_iterations == 0)
	    initial_L2_norm = L2_error;
	  
	  L2_error_relative = L2_error / initial_L2_norm;



	  if (bool_use_adaptive_newton_bound && (L2_error_relative > 1.0e-4))
	    {
	      lower_bound_newton_update = 1.0e-4 * L2_error_relative;
	      lower_bound_newton_residuum = 1.0e-4 * L2_error_relative; 
	    }
	  else
	    {
	      lower_bound_newton_update = 1.0e-8; 
	      lower_bound_newton_residuum = 1.0e-8; 

	      if (current_test_case.IsMiehe_Shear())
		{
		  lower_bound_newton_update = 1.0e-8;
		  lower_bound_newton_residuum = 1.0e-8; 
		}

	    }



	  std::cout << std::endl;
	  std::cout << "AL error: " 
		    << penal_iterations << "\t" 
		    << solution_lambda_penal_func.l2_norm() << "\t" 
		    << L2_error << "\t"
		    << L2_error_relative << "\t"
		    << lower_bound_newton_residuum
		    << std::endl;
	  std::cout << std::endl; 

	  
	  penal_iterations++; 

	}
      while ((L2_error > tolerance_absolute_augmented_L_iteration) &&
	     (L2_error_relative > tolerance_augmented_L_iteration) && 
	     (penal_iterations < max_no_of_augmented_L_penal_iterations));
      //      std::cout << std::endl;
      std::cout << "NumALIter: " << penal_iterations << std::endl;



  
      /*  
	  // TODO: evtl wieder einkommentieren
      // Option 1 (check at each time step first for the larger pressure increment)
      if ((test_case == "pressurized") && 
	  (penal_iterations == max_no_of_augmented_L_penal_iterations))
	{

	  solution = tmp_old_sol;
	  penal_iterations = 0;
	  current_pressure -= pressure_increment;
	  current_pressure = current_pressure + pressure_increment_factor * pressure_increment;
	  pressure_increment_factor *= 0.1; // TODO 0.1
	  if (pressure_increment_factor < 1e-6)
	    {
	      std::cout << "Aborting. Too small pressure_increment." << std::endl;
	      abort();
	    }
	  
	  std::cout << "Repeat step with smaller pressure increment: " << time << "\t" << current_pressure << "\t" << pressure_increment_factor << std::endl;

	  goto pf_extra_step;
	}
      */


/*
      // Option 2 (once smaller increment is used, keep this smaller one)
      if ((test_case == "pressurized") && 
	  (penal_iterations == max_no_of_augmented_L_penal_iterations))
	{

	  solution = tmp_old_sol;
	  penal_iterations = 0;
	  current_pressure -= pressure_increment;
	  pressure_increment =  pressure_increment_factor * pressure_increment
	  current_pressure = current_pressure + pressure_increment;

	  pressure_increment_factor *= 0.1;

	  if (pressure_increment_factor < 1e-6)
	    {
	      std::cout << "Aborting. Too small pressure_increment." << std::endl;
	      abort();
	    }
	  
	  std::cout << "Repeat step with smaller pressure increment: " << time << "\t" << pressure_increment << "\t" << pressure_increment_factor << std::endl;

	  goto pf_extra_step;
	}
*/            



      // Nov 21, 2016: Ist bei predictor-corrector refinement
      // und augmented Lagrangian scheinbar hier notwendig
      // Entscheidend war aber eine hoehere Anzahl von line search
      // steps (anstatt 4 nehme ich wieder 10)
      project_back_phase_field();
          
      // TODO
      // Predictor-correcter check
      if (current_test_case.IsMiehe_Shear() ||
	  current_test_case.IsL_Shaped())
	{
	  // Plot predictor solution on old mesh
	  if (bool_plot_additional_solutions)
	    {
	      output_results (prec_corr_counter_output_a,solution);
	      prec_corr_counter_output_a += 10000;
	    }

	  bool changed = refine_mesh();
	  if (changed)
	    {
	      // redo the current time step
	      std::cout << "Mesh change: corrector step" << std::endl;

	      // Plot predictor solution on new mesh
	      if (bool_plot_additional_solutions)
		{
		  output_results (prec_corr_counter_output_b,solution);
		  prec_corr_counter_output_b += 10000;
		}

	      //time -= timestep;
	      solution = old_timestep_solution;
	      
	      // Plot old solution on new mesh
	      if (bool_plot_additional_solutions)
		{
		  output_results (prec_corr_counter_output_c,solution);
		  prec_corr_counter_output_c += 10000;
		}

	      goto redo_step;

	    }
	}


	// Total number of Newton steps per time step
      std::cout << "NumNewtonIterTotal: " << total_number_newton_steps << std::endl;
      std::cout << std::endl;

}


