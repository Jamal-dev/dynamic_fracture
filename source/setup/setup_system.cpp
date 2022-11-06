#include "./../../include/dynamic_fracture.h"

// This function is similar to many deal.II tuturial steps.
template <int dim>
void Dynamic_Fracture_Problem<dim>::setup_system ()
{
  timer.enter_section("Setup system.");


  value_phase_field_for_refinement = 0.9;


  // We set runtime parameters to drive the problem.
  // These parameters could also be read from a parameter file that
  // can be handled by the ParameterHandler object (see step-19)
  system_matrix.clear ();
  
  dof_handler.distribute_dofs (fe);  
  DoFRenumbering::Cuthill_McKee (dof_handler);

  std::vector<unsigned int> block_component (dim+1+dim,0);
  block_component[dim] = 1;
  block_component[dim+1] = 2;
  block_component[dim+2] = 2;
  DoFRenumbering::component_wise (dof_handler, block_component);

  {				 
    constraints.clear ();
    set_newton_bc ();
    DoFTools::make_hanging_node_constraints (dof_handler,
					     constraints);
  }
  constraints.close ();
  
  std::vector<unsigned int> dofs_per_block (3);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);  
  const unsigned int n_u = dofs_per_block[0],
    n_c =  dofs_per_block[1], n_v = dofs_per_block[2];

  std::cout << "Cells:\t"
            << triangulation.n_active_cells()
            << std::endl  	  
            << "DoFs:\t"
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_c << '+' << n_v <<  ')'
            << std::endl;


 
      
 {
   BlockDynamicSparsityPattern csp (3,3);

    csp.block(0,0).reinit (n_u, n_u);
    csp.block(0,1).reinit (n_u, n_c);
    csp.block(0,2).reinit (n_u, n_v);
  
    csp.block(1,0).reinit (n_c, n_u);
    csp.block(1,1).reinit (n_c, n_c);
    csp.block(1,2).reinit (n_c, n_v);

    csp.block(2,0).reinit (n_v, n_u);
    csp.block(2,1).reinit (n_v, n_c);
    csp.block(2,2).reinit (n_v, n_v);

 
    csp.collect_sizes();    
  

    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);

    sparsity_pattern.copy_from (csp);
  }
 
 system_matrix.reinit (sparsity_pattern);

  // Actual solution at time step n
  solution.reinit (3);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_c);
  solution.block(2).reinit (n_v);
  
  solution.collect_sizes ();
 
  // Old timestep solution at time step n-1
  old_timestep_solution.reinit (3);
  old_timestep_solution.block(0).reinit (n_u);
  old_timestep_solution.block(1).reinit (n_c);
  old_timestep_solution.block(2).reinit (n_v);
 
  old_timestep_solution.collect_sizes ();

  // Old timestep solution at time step n-2
  old_old_timestep_solution.reinit (3);
  old_old_timestep_solution.block(0).reinit (n_u);
  old_old_timestep_solution.block(1).reinit (n_c);
  old_old_timestep_solution.block(2).reinit (n_v);
 
  old_old_timestep_solution.collect_sizes ();

  // temporary solution for time adaptivity 
  tmp_solution_for_time_adaptivity.reinit (3);
  tmp_solution_for_time_adaptivity.block(0).reinit (n_u);
  tmp_solution_for_time_adaptivity.block(1).reinit (n_c);
  tmp_solution_for_time_adaptivity.block(2).reinit (n_v);

  tmp_solution_for_time_adaptivity.collect_sizes ();


  // Updates for Newton's method
  newton_update.reinit (3);
  newton_update.block(0).reinit (n_u);
  newton_update.block(1).reinit (n_c);
  newton_update.block(2).reinit (n_v);

  newton_update.collect_sizes ();
 
  // Residual for  Newton's method
  system_rhs.reinit (3);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_c);
  system_rhs.block(2).reinit (n_v);

  system_rhs.collect_sizes ();


  // Lambda penal function
  solution_lambda_penal_func.reinit (3);
  solution_lambda_penal_func.block(0).reinit (n_u);
  solution_lambda_penal_func.block(1).reinit (n_c);
  solution_lambda_penal_func.block(2).reinit (n_v);

  solution_lambda_penal_func.collect_sizes ();

  old_timestep_solution_lambda_penal_func.reinit (3);
  old_timestep_solution_lambda_penal_func.block(0).reinit (n_u);
  old_timestep_solution_lambda_penal_func.block(1).reinit (n_c);
  old_timestep_solution_lambda_penal_func.block(2).reinit (n_v);

  old_timestep_solution_lambda_penal_func.collect_sizes ();


  // TODO: check, if this needs to be updated
  // after refining the mesh (currently it is because
  // setup_system is called in the refine_mesh function)
  if (bool_use_strain_history)
    setup_quadrature_point_history ();

  double min_cell_diameter = 1.0e+10;
 typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  for (; cell!=endc; ++cell)
    { 
      cell_diameter = cell->diameter();
      if (min_cell_diameter > cell_diameter)
	min_cell_diameter = cell_diameter;	

    }

  std::cout << "Min cell dia: " << min_cell_diameter << std::endl;

  timer.exit_section(); 
}


