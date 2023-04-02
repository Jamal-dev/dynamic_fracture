#include "./../../include/dynamic_fracture.h"



// In this method, we set up runtime parameters 
template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_p_mesh_1 ()
{
  // it does not have anything with boundary condition
  // changed bool_initial_crack_via_phase_field = false;
  bool_initial_crack_via_phase_field = false;
  

  // Parameters
  current_pressure = 0.0; 
  alpha_biot = 0.0;

  G_c = 1.5e3; // N/mm
  delta_penal = 0.0; // simple penalization
  gamma_penal = 1.0; //  augmented Lagrangian penalization

  gravity_x = 0.0;
  gravity_y = 0.0;
  volume_source = 0.0;


  // Elasticity parameters
  force_structure_x = 0.0;
  force_structure_y = 0.0;

  // Traction 
  traction_x = 0.0;
  traction_y = 0.0; 


  density_structure = 1.2e2; // kg/m^3 
  poisson_ratio_nu = 0.3;//0.3; 
  E_modulus = 2e11;//2e11; // pa
  // Timestep size:
  timestep = 1.0e-5;//1.0e-7; //1.0e-4;
  end_time_value = 0.00375;//1.0e-2; // Crack reaches lower left around 1.3e-2 sec
  // Structure parameters
  lame_coefficient_mu = E_modulus / (2.0 * (1 + poisson_ratio_nu));

  lame_coefficient_lambda = (2 * poisson_ratio_nu * lame_coefficient_mu)
                                / (1.0 - 2 * poisson_ratio_nu);
  

  
  // Timestepping schemes
  //BE, CN, CN_shifted
  // BE: backward euler
  time_stepping_scheme = "BE";

  

  // Maximum number of timesteps:
  max_no_timesteps = 1e6;//1e5; //130;
  // changed end_time_value = 1.0e-2;
  

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
  // Comented this part
  if (refinement_level>3)
    std::logic_error("Not implemented yet");
  std::string mesh_file_name = "mesh_files/example1/mesh_" + std::to_string(refinement_level) +".msh";
  if (!bool_initial_crack_via_phase_field)
    grid_name  = mesh_file_name;
  else if (bool_initial_crack_via_phase_field)
    std::logic_error("Not implemented yet"); 
 
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==2, ExcInternalError());
  grid_in.read_msh(input_file); 
  
  global_refinement_steps = 0;
  pred_corr_levels = 1;
  triangulation.refine_global (global_refinement_steps); 

  //filename_basis  = "solution_Miehe_eps_2h_ref_6_delta_0_"; 
  std::string parent_dir = "./results/patrick_tests/example1/mesh_" + std::to_string(refinement_level) ;
  parent_directory_ = parent_dir;
  create_directory (parent_dir);
  filename_basis  = parent_dir + "/" +  "solution_p_mesh" + std::to_string(refinement_level) +"_test_"; 
  bool_use_error_oriented_Newton = false;
  bool_use_modified_Newton = true; // if true need to set error_oriented_Newton to false
  bool_set_explicitely_delta_fp = false; // if true, must set use_modified Newton to false
  bool_set_explicitely_constant_k = false;
  bool_use_pf_extra = true;

  a_fp = 1e-2; //0.001; //0.05;
  b_fp = 5.0; //1.5;
  max_no_of_augmented_L_penal_iterations = 10; // 20
  if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

  // changed tolerance_augmented_L_iteration = 1.0e-5;
  tolerance_augmented_L_iteration = 1.0e-5;
  tolerance_absolute_augmented_L_iteration = 1.0e-5;

  max_no_newton_steps  = 100;

  // Line search
  max_no_line_search_steps = 10; 
  if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  // Redefined below !!
  lower_bound_newton_update = 1.0e-10;
  lower_bound_newton_residuum = 1.0e-10; 

  // When `true' be careful because tolerances are fixed 
  //  and differ from the choices here.
  bool_use_adaptive_newton_bound = false;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();

  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = true;

}


