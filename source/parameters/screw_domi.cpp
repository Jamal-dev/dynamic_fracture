#include "./../../include/dynamic_fracture.h"

// In this method, we set up runtime parameters 
template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_screw_domi ()
{
  // screw
  //sub_test_case = "hohlwalzung";
 


  // Parameters
  current_pressure = 0.0; 
  alpha_biot = 0.0;

  G_c = 2.7;
  delta_penal = 0.0; // simple penalization
  gamma_penal = 1.0; // augmented Lagrangian penalization

  gravity_x = 0.0;
  gravity_y = 0.0;
  volume_source = 0.0;


  // Elasticity parameters
  force_structure_x = 0.0;
  force_structure_y = 0.0;

  // Traction 
  traction_x = 0.0;
  traction_y = 0.0; 


  density_structure = 1.0; 

  // Solid parameters
  lame_coefficient_mu = 80.77e+3;
  lame_coefficient_lambda = 121.15e+3;


  
  // Timestepping schemes
  //BE, CN, CN_shifted
  time_stepping_scheme = "BE";

  // Timestep size:
  timestep = 1.0e-2; //1.0e-4;

  // Maximum number of timesteps:
  max_no_timesteps = 20; //130;
  end_time_value = 1.0; // Crack reaches lower left around 1.3e-2 sec

  number_of_nonadaptive_time_steps = 2;

  TOL_for_time_adaptivity = 1.0e-2; // 1.0e-2
  use_time_step_control   = false;


 
  timestep_rejection_factor = 0.5; // 0.1, 0.2 (Turek's suggestion)
  timestep_growth_factor = 1.0e+4; // 2.0
 
  timestep_upper_bound = 5.0e-5;
  timestep_lower_bound = 5.0e-6;



  
  // A variable to count the number of time steps
  timestep_number = 0;

  // Counts total time  
  time = 0;
 
  // Here, we choose a time-stepping scheme that
  // is based on finite differences:
  // BE         = backward Euler scheme 
  // CN         = Crank-Nicolson scheme
  // CN_shifted = time-shifted Crank-Nicolson scheme 
  // For further properties of these schemes,
  // we refer to standard literature.
  if (time_stepping_scheme == "BE")
    theta = 1.0;
  else if (time_stepping_scheme == "CN")
    theta = 0.5;
  else if (time_stepping_scheme == "CN_shifted")
    theta = 0.5 + timestep;
  else 
    std::cout << "No such timestepping scheme" << std::endl;

  // In the following, we read a *.inp grid from a file.
  std::string grid_name;
  grid_name  = "Schraube_2d.msh"; 
 
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==2, ExcInternalError());
  grid_in.read_msh (input_file); 

  typename Triangulation<dim>::cell_iterator
    cell = triangulation.begin(),
    endc = triangulation.end();
    for (; cell!=endc; ++cell)
      for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
	if (cell->at_boundary(i))
	  if (cell->face(i)->center()(1) == 0)
	    {
	      cell->face(i)->set_boundary_id(3);
	    }
	  else if (cell->face(i)->center()(1) == -17.40)
	    {
	      cell->face(i)->set_boundary_id(2);
	    }
    
  
  global_refinement_steps = 1;
  pred_corr_levels = 0;
  triangulation.refine_global (global_refinement_steps); 

  filename_basis  = "solution_screw_domi_"; 
  //filename_basis  = "solution_test_"; 
  bool_use_error_oriented_Newton = false;
  bool_use_modified_Newton = true; // if true and be used, then must set error oriented false
  bool_set_explicitely_delta_fp = false; // if true, must set use_modified Newton to false
  bool_set_explicitely_constant_k = false;
  bool_use_pf_extra = false;

  a_fp = 0.01; // 0.01 for modified Newton
  b_fp = 5.0; // 5.0 for modified Newton
  max_no_of_augmented_L_penal_iterations = 50; // 50
  if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

 tolerance_augmented_L_iteration = 1.0e-3; //1.0e-4;
  tolerance_absolute_augmented_L_iteration = 1.0e-3;

  max_no_newton_steps  = 100;
  max_no_line_search_steps = 10;
 if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  lower_bound_newton_update = 1.0e-8;
  lower_bound_newton_residuum = 1.0e-8; 

  // When `true' be careful because tolerances are fixed 
  //  and differ from the choices here.
  bool_use_adaptive_newton_bound = false;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();

  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = false;

}

