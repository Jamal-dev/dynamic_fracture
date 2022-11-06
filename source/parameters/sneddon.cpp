#include "./../../include/dynamic_fracture.h"


// In this method, we set up runtime parameters 
template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_Sneddon ()
{

  // Parameters
  current_pressure = 1.0e-3; 
  alpha_biot = 0.0;

  G_c = 1.0;
  delta_penal = 0.0; // simple penalization
  gamma_penal = 1.0e+3; //1.0e+3; //1.0e+4; // augmented Lagrangian penalization

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
  density_structure = 1.0; 
  poisson_ratio_nu = 0.2; 
  E_modulus = 1.0;
  E_prime = E_modulus/(1.0 - poisson_ratio_nu * poisson_ratio_nu);
  
  lame_coefficient_mu = E_modulus/(2.0*(1 + poisson_ratio_nu));
  
  lame_coefficient_lambda = (2 * poisson_ratio_nu * lame_coefficient_mu)/(1.0 - 2 * poisson_ratio_nu);

  
  // Timestepping schemes
  //BE, CN, CN_shifted
  time_stepping_scheme = "BE";

  // Timestep size:
  timestep = 1.0; //2.5e-1; //1.0;

  // Maximum number of timesteps:
  max_no_timesteps = 5;
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
  //grid_name  = "unit_slit.inp"; 
  // Example 2
  grid_name  = "unit_square_4.inp"; 
  
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==2, ExcInternalError());
  grid_in.read_ucd (input_file); 


 // TODO: insert all time step tolerances etc.
  
  global_refinement_steps = 6;
  pred_corr_levels = 2;
  triangulation.refine_global (global_refinement_steps); 

  //filename_basis  = "solution_pressurized_pf_extra_eps_2h_ref_7_adaptive_"; 
  filename_basis  = "solution_Sneddon_pf_extra_2h_ref_5_"; 
  filename_basis_cod = "cod_ref_6_";
  bool_use_error_oriented_Newton = false;
  bool_use_pf_extra = true;
  max_no_of_augmented_L_penal_iterations = 50; //15; //40;
 if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

  tolerance_augmented_L_iteration = 1.0e-3; // 1e-4;
  tolerance_absolute_augmented_L_iteration = 1.0e-8;

  max_no_newton_steps  = 20; //20; // 50
  max_no_line_search_steps = 10; 
  if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  lower_bound_newton_update = 1.0e-8;
  lower_bound_newton_residuum = 1.0e-8; 

  bool_use_adaptive_newton_bound = true;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();

  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = false;


}


