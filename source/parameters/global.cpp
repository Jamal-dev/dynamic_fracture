#include "./../../include/dynamic_fracture.h"

template <int dim>
void Dynamic_Fracture_Problem<dim>::set_global_parameters ()
{
  // Initialization of some variables
  // that may be overwritten in the set_runtime routines
  bool_set_explicitely_delta_fp = false;
  bool_set_explicitely_constant_k = false;
  
  bool_use_strain_history = false;
  bool_initial_crack_via_phase_field = false;

  upper_newton_rho = 0.9;

  constant_k = 1.1e-2; // it's kappa
  alpha_eps = 2.2e-2; // it's eps

  use_time_step_control   = false;

  bool_use_dynamic_code_with_velocities = true;

}
