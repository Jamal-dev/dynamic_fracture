#ifndef dynamic_fracture_local_h
#define dynamic_fracture_local_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>  

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>

// #include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

//-----------------------------------------------------------//
#include <sys/types.h>
#include <dirent.h>
#include <cstdlib>

#include "test_cases.h"

void create_directory (const std::string &directory_name)
{
  DIR *dir;
  if ((dir = opendir(directory_name.c_str())) == NULL)
    {
      std::string command = "mkdir -p " + directory_name;
      system(command.c_str());
      std::cout<<"Created directory "<<directory_name<<std::endl;
    }
  else
    closedir(dir);
}

//------------------------------------------------------------//
 template <int dim>
  struct PointHistory
  {
    Tensor<2,dim> old_stress;
  };


/////////////////////////////////////////////////////////////////////

template <int dim>
class Dynamic_Fracture_Problem 
{
public:
  
  Dynamic_Fracture_Problem (const unsigned int degree);
  ~Dynamic_Fracture_Problem (); 
  void run ();
  
//private:
  
  void set_global_parameters ();

  void set_runtime_parameters_p_mesh1 ();
  void set_runtime_parameters_Miehe ();
  void set_runtime_parameters_Dynamic_Slit ();
  void set_runtime_parameters_L_shaped ();
  void set_runtime_parameters_Sneddon ();
  void set_runtime_parameters_pressurized ();
  void set_runtime_parameters_screw_domi ();
  void set_runtime_parameters_Sneddon3D ();
  void set_runtime_parameters_Het3D ();

  void setup_system ();
  void make_material_vectors();
  void assemble_system_matrix ();   
  void assemble_system_rhs ();
  
  // functions for setting initial and boundary condition
  void set_initial_bc (const double time);
  void set_newton_bc ();
  
  void solve ();
  double newton_iteration(const double time);

  void solve_spatial_problem();

 // Nonlinear solver: error-based descend
  void newton_iteration_error_based(const double time);	

		  
  void output_results (const unsigned int refinement_cycle,
		       const BlockVector<double> solution) const;
  
  // Compute functional values
  double compute_point_value (Point<dim> p,
			      const unsigned int component) const;
  
  void compute_functional_values ();
 

  void compute_functional_values_Sneddon (); 
  double compute_cod (const double eval_line); 
  void compute_cod_Sneddon3D();

  void compute_stress_per_cell ();

  void compute_energy ();

  double goal_functional_stress_x ();

  // Mesh refinement
  bool refine_mesh();
  void project_back_phase_field (); 


  void setup_quadrature_point_history ();
  void update_quadrature_point_history ();

  std::vector<PointHistory<dim> >  quadrature_point_history;

  const unsigned int   degree;
  
  Triangulation<dim>   triangulation;
  FESystem<dim>        fe;
  DoFHandler<dim>      dof_handler;

  ConstraintMatrix     constraints;
  
  BlockSparsityPattern      sparsity_pattern; 
  BlockSparseMatrix<double> system_matrix; 
  
  BlockVector<double> solution, newton_update, old_timestep_solution, 
    old_old_timestep_solution, tmp_solution_for_time_adaptivity;
  BlockVector<double> system_rhs;
  BlockVector<double> solution_lambda_penal_func, old_timestep_solution_lambda_penal_func;


  Vector<double> lame_coefficient_mu_vector, lame_coefficient_lambda_vector;
  Vector<float> solution_stress_per_cell;

  TimerOutput         timer;
  // ConditionalOStream pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_com) == 0));

  // MPI_Comm mpi_com;

  Function<dim> *func_emodulus;
  
  // Global variables for timestepping scheme   
  unsigned int timestep_number;
  unsigned int max_no_timesteps;  
  double timestep, theta, time; 
  std::string time_stepping_scheme;
  std::string test_case, sub_test_case;
  test_cases current_test_case;
  double old_timestep, old_old_timestep;

  double force_structure_x, force_structure_y;	  
  

  double gravity_x, gravity_y, volume_source, traction_x, traction_y;
  
  
  // Structure parameters
  double density_structure; 
  double density_timestep_ratio;
  double lame_coefficient_mu, lame_coefficient_lambda, poisson_ratio_nu;  

  // test cases for functional values function
  std::set<test_cases> remaining_test_cases_functional_values;
  double cell_diameter;  
  double max_no_of_augmented_L_penal_iterations;
  unsigned int penal_iterations;
  
  bool bool_use_stress_splitting;
 
  double constant_k, alpha_eps, G_c, delta_penal, gamma_penal;
  double upper_newton_rho;
  double tolerance_augmented_L_iteration,  tolerance_absolute_augmented_L_iteration;
  double current_pressure, old_timestep_current_pressure, alpha_biot;
  double E_modulus, E_prime;
  unsigned int total_number_newton_steps;


  double lower_bound_newton_update;
  double lower_bound_newton_residuum; 

  bool bool_use_pf_extra, bool_use_error_oriented_Newton, bool_use_modified_Newton, 
    bool_use_adaptive_newton_bound, bool_plot_additional_solutions,
    bool_set_explicitely_delta_fp, bool_set_explicitely_constant_k;
  unsigned int max_no_newton_steps, max_no_line_search_steps;
  std::string filename_basis, filename_basis_cod;

  double value_phase_field_for_refinement;

  unsigned int global_refinement_steps;
  unsigned int pred_corr_levels;

  double TOL_for_time_adaptivity, timestep_rejection_factor, timestep_growth_factor;
  bool use_time_step_control;
  double timestep_upper_bound, timestep_lower_bound,end_time_value;
  int number_of_nonadaptive_time_steps;

  double delta_fixed_point_newton;
  double a_fp, b_fp;

  bool bool_set_initial_strain_history, bool_initial_crack_via_phase_field, bool_use_strain_history;

  bool bool_use_dynamic_code_with_velocities;
  

};


/*
        Constructor and desturctor part
*/
// The constructor of this class is comparable 
// to other tutorials steps, e.g., step-22, and step-31.
template <int dim>
Dynamic_Fracture_Problem<dim>::Dynamic_Fracture_Problem (const unsigned int degree)
                :
                degree (degree),
		triangulation (Triangulation<dim>::maximum_smoothing),
                /*dim is the number of copies*/fe (FE_Q<dim>(degree), dim,  // velocities                
		    FE_Q<dim>(degree), 1,    // phase-field
		    FE_Q<dim>(degree), dim),   // displacements
                dof_handler (triangulation),
    timer (std::cout, TimerOutput::summary, TimerOutput::cpu_times)
		// timer (std::cout, TimerOutput::summary, TimerOutput::cpu_times)
    // timer (pcout, TimerOutput::summary, TimerOutput::cpu_and_wall_times)
    // timer(mpi_com, pcout, TimerOutput::every_call_and_summary,
    //     TimerOutput::cpu_and_wall_times)		
{}


// This is the standard destructor.
template <int dim>
Dynamic_Fracture_Problem<dim>::~Dynamic_Fracture_Problem () 
{}

#endif
