#include "./../../include/dynamic_fracture.h"

template <int dim>
void Dynamic_Fracture_Problem<dim>::newton_iteration_error_based (const double time) 
					       
{ 
   Timer timer_newton;
   //const double lower_bound_newton_update = 1.0e-8; //1.0e-12; 
   //const double lower_bound_newton_residuum = 1.0e-8; 
   //const unsigned int max_no_newton_steps  = 50;

  const double lambda_newton_lower_bound = 1.0e-10;

  // Decision whether the system matrix should be build
  // at each Newton step
  const double nonlinear_rho = 0.1; 
 
  // Line search parameters
  unsigned int line_search_step;
  const unsigned int  max_no_line_search_steps = 10;
  const double line_search_damping = 0.6;
  double new_newton_residuum;

  if (bool_set_explicitely_delta_fp)
    delta_fixed_point_newton = 0.0;

  
  // Application of the initial boundary conditions to the 
  // variational equations:
  //set_initial_bc (time);
  assemble_system_rhs();

  double newton_residuum = system_rhs.l2_norm(); 
  double newton_update_norm = 1.0;
  double old_newton_update_norm = newton_update_norm;
  double bar_newton_update_norm = newton_update_norm;
  double old_bar_newton_update_norm = newton_update_norm;
  unsigned int newton_step = 1;

  BlockVector<double> bar_newton_update = newton_update;
  BlockVector<double> old_bar_newton_update = bar_newton_update;
  BlockVector<double> new_newton_update = newton_update;
  double theta_newton = 1.0;
  double old_theta_newton = theta_newton;

  double lambda_newton = 1.0;
  double lambda_newton_prime = lambda_newton;
  double old_lambda_newton = lambda_newton;
  double tmp_lambda_newton = lambda_newton;

  double mu_newton = 1.0;
  double mu_newton_prime = 1.0;

  while (newton_update_norm > lower_bound_newton_update &&
	 newton_residuum > lower_bound_newton_residuum &&
	 newton_step < max_no_newton_steps)
    {
      timer_newton.start();

      // Step 1
      old_newton_update_norm = newton_update_norm;
      old_theta_newton = theta_newton;

      old_bar_newton_update = bar_newton_update;
      old_bar_newton_update_norm = bar_newton_update_norm;

      old_lambda_newton = lambda_newton;

      
      assemble_system_rhs();
 
      newton_residuum = system_rhs.l2_norm();

      if (newton_residuum < lower_bound_newton_residuum)
	{
	  std::cout << '\t' 
		    << std::scientific 
		    << newton_residuum << std::endl;
	  break;
	}





      //if (theta_newton > nonlinear_rho)
	assemble_system_matrix ();	


      // Solve Ax = b
      solve ();	

      new_newton_update = newton_update;
      newton_update_norm = new_newton_update.l2_norm(); 
      //      std::cout << newton_update_norm << std::endl;

      if (newton_update_norm < lower_bound_newton_update)
	{
	  solution += new_newton_update;
	  std::cout << '\t' 
		    << std::scientific 
		    << newton_update_norm << std::endl;
	  break;
	}

      

      BlockVector<double> diff_newton = bar_newton_update;
      for (unsigned int i=0; i < new_newton_update.size();i++)
	diff_newton[i] = old_bar_newton_update[i] - new_newton_update[i];

      double diff_newton_norm = diff_newton.l2_norm();

      mu_newton = 1.0;
      if (newton_step > 1)
      	mu_newton = (old_newton_update_norm * old_bar_newton_update_norm) / 
      	  (diff_newton_norm * newton_update_norm) * old_lambda_newton;
      

      lambda_newton = std::min(1.0, mu_newton);


    regu_test:
      //std::cout << "lambda newton: " << lambda_newton << std::endl;
      if (lambda_newton < lambda_newton_lower_bound)
	{
	  std::cout << lambda_newton << "\tConvergence failure. Aborting in Step 1." << std::endl;

	  // Option 1
	  //abort();	  

	  /*
	  // Option 2: switch to residual-based Newton
	  std::cout << lambda_newton << "\tConvergence failure in error-based Newton. Go to residual-based Newton." << std::endl;

	  // TODO: be careful, in standard newton the 
	  // initial bc are commented!!!!!
	  newton_iteration (time);
	  break;
	  */	  

	  
	  // Option 3: just take the small step and continue
	  std::cout << lambda_newton << "\tConvergence failure in error-based Newton. Take the smallest update and continue." << std::endl;
	  for (unsigned int i=0; i < new_newton_update.size();i++)
	    {
	      solution(i) = solution(i) + lambda_newton * new_newton_update(i);
	    }
	  break;	
	  

	}

      // Step 2
    step_2:

      tmp_lambda_newton = lambda_newton;
      for (unsigned int i=0; i < new_newton_update.size();i++)
	{
	  solution(i) = solution(i) + lambda_newton * new_newton_update(i);
	}

      assemble_system_rhs ();

      // Solve now the bar-system
      // old Jacobian and new right hand side
      solve ();
      bar_newton_update = newton_update;
      bar_newton_update_norm = bar_newton_update.l2_norm(); 
   
      // Step 3
      // Compute monitoring quantities
      BlockVector<double> diff_mu_newton_prime = new_newton_update;
      for (unsigned int i=0; i < new_newton_update.size();i++)
	{
	  diff_mu_newton_prime(i) =
	    bar_newton_update(i) - (1.0 - lambda_newton) * new_newton_update(i);
	}


      double diff_mu_newton_prime_norm = diff_mu_newton_prime.l2_norm();
 
      theta_newton = bar_newton_update_norm / newton_update_norm;
      //std::cout << "theta_newton: " <<  theta_newton << std::endl;


      mu_newton_prime = 1.0;
      if (newton_step > 1)
	mu_newton_prime = 0.5 * (newton_update_norm * lambda_newton * lambda_newton) /
	  diff_mu_newton_prime_norm;

      lambda_newton_prime = lambda_newton;
      double eps_theta_newton = 1.0; // 1.001
      if (theta_newton >= eps_theta_newton) //TODO 1.0
       {
	 lambda_newton_prime = std::min(mu_newton_prime, 0.5 * lambda_newton);
	 lambda_newton = lambda_newton_prime;
	 //std::cout << "Monitoring test: " << lambda_newton << std::endl;

	 for (unsigned int i=0; i < new_newton_update.size();i++)
	   {
	     solution(i) = solution(i) - tmp_lambda_newton * new_newton_update(i);
	   }
	 

	 goto regu_test;
	 
       }
     else 
       lambda_newton_prime = std::min(mu_newton_prime, 1.0);
        

      
     if ((lambda_newton < (lambda_newton_prime + 1.0e-10)) && 
	 (lambda_newton > (lambda_newton_prime - 1.0e-10)) && 
	 (lambda_newton < 1.0 + 1.0e-10) &&
	 (lambda_newton > 1.0 - 1.0e-10)
	 )
       {
//	 std::cout << "Choice 1: " 
//		   << lambda_newton << "  " 
//		   << bar_newton_update_norm << "  "
//		   << theta_newton << std::endl;

	 if (bar_newton_update_norm < lower_bound_newton_update)
	   {
	     solution += bar_newton_update;
	     std::cout << '\t' 
		       << std::scientific 
		       << bar_newton_update_norm << std::endl;
	     break;
	   }
	 
       }
//     else if (lambda_newton_prime >= 4.0 * lambda_newton)
//       {
//	 // std::cout << "Choice 2:" << lambda_newton_prime << "   "  << lambda_newton  << std::endl;
//
//	 for (unsigned int i=0; i < new_newton_update.size();i++)
//	   {
//	     solution(i) = solution(i) - tmp_lambda_newton * new_newton_update(i);
//	   }
//	 
//	 
//	 lambda_newton = lambda_newton_prime;
//
//	 goto step_2;
//       }
        
     
     //std::cout << "Choice 3: " << lambda_newton << std::endl;

      assemble_system_rhs();
      newton_residuum = system_rhs.l2_norm();

      timer_newton.stop();
      
      std::cout << std::setprecision(5) <<newton_step << '\t' 
		<< std::scientific << newton_update_norm << '\t'
		<< std::scientific << bar_newton_update_norm << '\t'
		<< std::scientific << newton_residuum << '\t'
		<< std::scientific << newton_update_norm/old_newton_update_norm  <<'\t' ;
      if (old_theta_newton > nonlinear_rho)
	std::cout << "r" << '\t' ;
      else 
	std::cout << " " << '\t' ;
        
      std::cout << delta_fixed_point_newton  << '\t'
		<< std::scientific << timer_newton.wall_time ()
		<< std::endl;


	if (bool_use_modified_Newton)
	{
	  // Update delta for dynamic switch between fixed point and Newton
	  double Qn = bar_newton_update_norm/newton_update_norm;
	  double Qn_inv = newton_update_norm/bar_newton_update_norm;
	  
	  // Mandel's et al. formula (31)
	  //delta_fixed_point_newton = delta_fixed_point_newton * (0.2 + 4.0/(0.7 + std::exp(1.5 * Qn)));
	  
	  delta_fixed_point_newton = delta_fixed_point_newton * (a_fp/(std::exp(Qn_inv)) + b_fp/(std::exp(Qn)));

	  // Normalize delta
	  if (delta_fixed_point_newton > 1.0)
	    delta_fixed_point_newton = 1.0;
	  else if (delta_fixed_point_newton < 0.0)
	    delta_fixed_point_newton = 0.0;
	  
	}




      // Updates
      timer_newton.reset();
      newton_step++;      
    }

  //std::cout << "NumNewtonIter: " << newton_step << std::endl;
  total_number_newton_steps += newton_step;

}
