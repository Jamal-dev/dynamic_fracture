#include "./../../include/dynamic_fracture.h"

// In this method, we set up runtime parameters 
template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_L_shaped ()
{

  // Parameters
  current_pressure = 0.0; 
  alpha_biot = 0.0;

  G_c = 8.9e-2;
  delta_penal = 0.0; // simple penalization
  gamma_penal = 1.0e-3; // augmented Lagrangian penalization

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

  // Structure parameters
  lame_coefficient_mu = 10.95e+3; 
  poisson_ratio_nu = 0.3; 
  
  lame_coefficient_lambda = 6.16e+3; //(2 * poisson_ratio_nu * lame_coefficient_mu)/(1.0 - 2 * poisson_ratio_nu);

  
  // Timestepping schemes
  //BE, CN, CN_shifted
  time_stepping_scheme = "BE";

  // Timestep size:
  timestep = 1.0e-3;

  // Maximum number of timesteps:
  max_no_timesteps = 2000;
  end_time_value = 1e+10;
 
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
  grid_name  = "l_shape_Jan_7_2016.inp"; 
  //grid_name  = "l_shape_Jan_8_2016.inp"; 
 
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==2, ExcInternalError());
  grid_in.read_ucd (input_file); 
  
  // TODO: insert all time step tolerances etc from Miehe oben
  // TODO: adjust goal functional
 
  global_refinement_steps = 2;
  pred_corr_levels = 0;   
  triangulation.refine_global (global_refinement_steps); 

  filename_basis  = "solution_l_shape_"; 
  bool_use_error_oriented_Newton = true;
  bool_use_pf_extra = true;

  max_no_of_augmented_L_penal_iterations = 20; // 50
 if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

  tolerance_augmented_L_iteration = 1.0e-4;
  tolerance_absolute_augmented_L_iteration = 1.0e-8;

  max_no_newton_steps  = 20; //50;
  max_no_line_search_steps = 10;
  if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  // reset in run method because of adaptive Newton technique
  lower_bound_newton_update = 1.0e-8;
  lower_bound_newton_residuum = 1.0e-8; 

  bool_use_adaptive_newton_bound = false;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();
 
  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = true;

}


