#include "./../../include/dynamic_fracture.h"

// This is the Newton iteration to solve the 
// non-linear system of equations. First, we declare some
// standard parameters of the solution method. Addionally,
// we also implement an easy line search algorithm. 
template <int dim>
double Dynamic_Fracture_Problem<dim>::newton_iteration (const double time) 
					       
{ 
  Timer timer_newton;
  //const double lower_bound_newton_residuum = 1.0e-8; 
  //const unsigned int max_no_newton_steps  = 20;

  // Decision whether the system matrix should be build
  // at each Newton step
  const double nonlinear_rho = 0.1; 
 
  // Line search parameters
  unsigned int line_search_step;

  // For residual-based Newton with increasing
  // residual use a small number of steps: 4
  //const unsigned int  max_no_line_search_steps = 10;
  const double line_search_damping = 0.6;
  double new_newton_residuum;
  
  // Application of the initial boundary conditions to the 
  // variational equations:
  // Removed because of coupled Newton iteration.
  assemble_system_rhs();

  double newton_residuum = system_rhs.linfty_norm(); 
  double old_newton_residuum= newton_residuum;
  unsigned int newton_step = 1;

  if (bool_set_explicitely_delta_fp)
    delta_fixed_point_newton = 0.0;
   
  if (newton_residuum < lower_bound_newton_residuum)
    {
      std::cout << '\t' 
		<< std::scientific 
		<< newton_residuum 
		<< std::endl;     
    }
  
  while (newton_residuum > lower_bound_newton_residuum &&
	 newton_step < max_no_newton_steps)
    {
      timer_newton.start();

  if (bool_set_explicitely_constant_k)
    {
      //delta_fixed_point_newton = 0.0;
      //delta_fixed_point_newton = newton_residuum;
      constant_k = newton_residuum;
      if (constant_k > (0.5 * alpha_eps))
	constant_k = 0.5 * alpha_eps;
      //std::cout << constant_k << std::endl;
    }
 

      old_newton_residuum = newton_residuum;
      
      assemble_system_rhs();
      newton_residuum = system_rhs.linfty_norm();

      if (newton_residuum < lower_bound_newton_residuum)
	{
	  std::cout << '\t' 
		    << std::scientific 
		    << newton_residuum << std::endl;
	  break;
	}

      if (newton_residuum > 1e+14)
	{
	  std::cout << '\t' 
		    << std::scientific 
		    << newton_residuum << std::endl;
	  std::cout << "Newton residual too high. Aborting." << std::endl;
	  abort();
	}


      // TODO: This seems to lead to 
      // non-robust Newton iterations
      //if (newton_residuum/old_newton_residuum > nonlinear_rho)
	assemble_system_matrix ();	

      // Solve Ax = b
      solve ();	  
        
      line_search_step = 0;	  
      for ( ; 
	    line_search_step < max_no_line_search_steps; 
	    ++line_search_step)
	{	 
	  solution +=newton_update;

	  assemble_system_rhs ();			
	  new_newton_residuum = system_rhs.linfty_norm();
	  
	  if (new_newton_residuum < newton_residuum)
	      break;
	  else 	  
	    solution -= newton_update;
	  
	  newton_update *= line_search_damping;

	}

      // Allow for increasing residual
      if (line_search_step == max_no_line_search_steps)
	{
	  solution +=newton_update;	   

	}
     
      timer_newton.stop();
      
      std::cout << std::setprecision(5) <<newton_step << '\t' 
		<< std::scientific << newton_residuum << '\t'
		<< std::scientific << newton_residuum/old_newton_residuum  <<'\t' ;
      if (newton_residuum/old_newton_residuum > nonlinear_rho)
	std::cout << "r" << '\t' ;
      else 
	std::cout << " " << '\t' ;
      std::cout << line_search_step  << '\t' 
		<< delta_fixed_point_newton  << '\t'
		<< std::scientific << timer_newton.wall_time ()
		<< std::endl;

//      if ((newton_residuum/old_newton_residuum > upper_newton_rho) && (newton_step > 1))
//	{
//	  break;
//	}


      if (bool_use_modified_Newton)
	{
	  // Update delta for dynamic switch between fixed point and Newton
	  double Qn = newton_residuum/old_newton_residuum;
	  double Qn_inv = old_newton_residuum/newton_residuum;
	  
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

  return newton_residuum/old_newton_residuum; 
}

