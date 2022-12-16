// Phase-field fracturetimer_newton
// Properties:
// - computes fracture propagation in elasticity
// - fully-coupled displacement phase-field model
// - predictor-corrector mesh adaptivity for small eps
// - tuncation-based goal-oriented adaptive time step control

// TODO: goal functionals
// augmented Lagrangian stopping only phase-field




#include "./include/utils.h"
// Importing initial values
#include "include/initial_values/pressurized.h"
#include "include/initial_values/sneddon3d.h"
#include "include/initial_values/het3d.h"
#include "include/initial_values/miehe.h"
#include "include/initial_values/hohlwalzung.h"
#include "include/initial_values/dynamic_slit.h"

// Importing boundary values
#include "source/boundary_conditions/nonhomoDirichletBoundary.cpp"
#include "source/boundary_conditions/nonhomoDirichletBoundaryVelocities.cpp"

// Importing math operations
#include "source/operations.cpp"

// Importing dynamic fracture model
#include "include/dynamic_fracture.h"

// Import runtime parameters
#include "source/parameters/global.cpp" 
#include "source/parameters/miehe.cpp"
#include "source/parameters/pressurized.cpp"
#include "source/parameters/sneddon3d.cpp"
#include "source/parameters/sneddon.cpp"
#include "source/parameters/het3d.cpp"
#include "source/parameters/screw_domi.cpp"
#include "source/parameters/l_shaped.cpp"
#include "source/parameters/dynamic_slit.cpp"
#include "source/parameters/p_mesh_1.cpp"
#include "source/parameters/p_notched_cavity.cpp"


// Import initial guesses for the newton iteration
#include "source/boundary_conditions/initial_guess_newton.cpp"
#include "source/boundary_conditions/newton_bc.cpp"

// Import assemble
#include "source/assemble/assemble_lhs.cpp"
#include "source/assemble/assemble_rhs.cpp"

// Import setupsystem
#include "source/setup/setup_system.cpp"
#include "source/setup/setup_quadrature_point_history.cpp"

// Import solve
#include "source/solve/solve.cpp"
#include "source/solve/back_phase_field.cpp"
#include "source/solve/newton_iteration.cpp"
#include "source/solve/newton_iteration_error_based.cpp"
#include "source/solve/solve_spatial_problem.cpp"

// Import output
#include "source/output/write_vtk.cpp"

// Import compute quanties of interest
#include "source/compute_quantities/cod.cpp"
#include "source/compute_quantities/cod_sneddon3d.cpp"
#include "source/compute_quantities/compute_energy.cpp"
#include "source/compute_quantities/functional_values.cpp"
#include "source/compute_quantities/functional_values_sneddon.cpp"
#include "source/compute_quantities/goal_functional_stress.cpp"
#include "source/compute_quantities/point_values.cpp"
#include "source/compute_quantities/stress_per_cell.cpp"

// Import mesh adaptivity
#include "source/mesh/refine_mesh.cpp"

// Update quadrature point history
#include "source/update_quantities/update_quadrature_point_history.cpp"

// Create quantities
#include "source/create_quantities/make_material_vectors.cpp"

// At the end of this top-matter, we import
// all deal.II names into the global
// namespace:				
using namespace dealii;
using namespace std;



////////////////////////////////////////////////////////////


template<int dim> void Dynamic_Fracture_Problem<dim>::write_rutime_parameters_csv()
{
   try
        {
            std::string filename = parent_directory_+  "/"+"parameters_current_simulation.csv";
            csvfile csv(filename,","); // throws exceptions!
            // Header
            csv.AddRow("Cells",triangulation.n_active_cells(),"DOFs",dof_handler.n_dofs()); 
            csv.AddRow("Dofs_u",get<0>(dofs_per_component) ,"Dofs_p",get<1>(dofs_per_component) ,"Dofs_v",get<2>(dofs_per_component) );
            csv.AddRow("lambda",lame_coefficient_lambda,"mu",lame_coefficient_mu);
            csv.AddRow("density",density_structure, "G_c", G_c);
            csv.AddRow("E",E_modulus,"nu",poisson_ratio_nu);
            csv.AddRow("timestep",timestep,"end_time",end_time_value);
            csv.AddRow("kappa",constant_k,"eps",alpha_eps);
            csv.AddRow("min_h", min_cell_diameter, "delta_penal",delta_penal);
            csv.AddRow("is_use_stress_splitting", bool_use_stress_splitting, "is_use_adaptive_newton_bound",bool_use_adaptive_newton_bound);
            
        }
    catch (const std::exception& ex)
        {std::cout << "Exception was thrown: " << ex.what() << std::endl;}  

}

// As usual, we have to call the run method. It handles
// the output stream to the terminal.
// Second, we define some output skip that is necessary 
// (and really useful) to avoid to much printing 
// of solutions. For large time dependent problems it is 
// sufficient to print only each tenth solution. 
// Third, we perform the time stepping scheme of 
// the solution process.
template <int dim>
void Dynamic_Fracture_Problem<dim>::run () 
{ 
  // Switch dimension !!
  current_test_case = test_cases::P_MESH_1;
  refinement_level = 1;
  // current_test_case = test_cases::P_NOTCHED_CAVITY;
  // Defining test cases
  // test_case = "dynamic_slit";
  // test_case = "miehe_shear";
  // test_case = "miehe_tension";
  //test_case = "l_shaped";
  //test_case = "Sneddon";
  // before it was pressurized test case
  // test_case = "pressurized";
  //test_case = "screw_domi";
  //test_case = "Sneddon3D";
  //test_case = "Het3D";

  // Switch dimension !!
  

  // Initialization of some variables
  // that may be overwritten in the set_runtime routines
  set_global_parameters ();

  // Setting specific parameters for each problem
  if (current_test_case == test_cases::MIEHE_TENSION ||
      current_test_case == test_cases::MIEHE_SHEAR)
    set_runtime_parameters_Miehe ();
  else if (current_test_case == test_cases::DYNAMIC_SLIT)
    set_runtime_parameters_Dynamic_Slit();
  else if (current_test_case == test_cases::P_MESH_1)
    set_runtime_parameters_p_mesh_1();
	else if (current_test_case == test_cases::P_NOTCHED_CAVITY)
    set_runtime_parameters_p_notched_cavity();
  else if (current_test_case == test_cases::L_SHAPED)
    set_runtime_parameters_L_shaped ();
  else if (current_test_case == test_cases::SNEDDON)
    set_runtime_parameters_Sneddon ();
  else if (current_test_case == test_cases::PRESSURIZED)
    set_runtime_parameters_pressurized ();
  else if (current_test_case == test_cases::SCREW_DOMI)
    set_runtime_parameters_screw_domi ();
  else if (current_test_case == test_cases::SNEDDON3D)
    {
      set_runtime_parameters_Sneddon3D ();
      // TODO: set dimension to 3 !!
    }
  else if (current_test_case == test_cases::HET3D)
    {
      set_runtime_parameters_Het3D ();
      // TODO: set dimension to 3 !!
    }
  else 
    {
      std::cout << "Framework not implemented. Aborting. " << std::endl;
      abort();
    }
  
  setup_system();

  write_rutime_parameters_csv();


  if (current_test_case == test_cases::SNEDDON3D 
      || current_test_case == test_cases::HET3D
      )
    {
      for (unsigned int i=0;i<3;i++)
	        refine_mesh();
    }


  min_cell_diameter = 1.0e+10;
  /*
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  for (; cell!=endc; ++cell)
    { 
      cell_diameter = cell->diameter();
      if (min_cell_diameter > cell_diameter)
	min_cell_diameter = cell_diameter;	
      
    }
  */
  
  min_cell_diameter = dof_handler.begin(0)->diameter();
  min_cell_diameter *= std::pow(2.,-1.0*( global_refinement_steps  +  pred_corr_levels ));
  std::cout << "Min cell dia: " << min_cell_diameter << std::endl;


  // Set phase-field parameters
  constant_k = 1.0e-12;// * min_cell_diameter;
  // TODO - adjust according to mesh levels
  alpha_eps =  2.0 * min_cell_diameter; 

  //if (test_case == "miehe_shear")
  //  alpha_eps =  0.0884;
  

  std::cout << "\n==============================" 
	    << "====================================="  << std::endl;
   std::cout << "Parameters\n" 
	    << "==========\n"
	     << "min h:             "   << min_cell_diameter << "\n"
	     << "k:                 "   <<  constant_k << "\n"  
	     << "eps:               "   <<  alpha_eps << "\n"
	     << "G_c:               "   <<  G_c << "\n"
	     << "delta penal:       "   <<  delta_penal << "\n"
	     << "Poisson nu:        "   <<  poisson_ratio_nu << "\n"
	     << "Lame mu:           "   <<  lame_coefficient_mu << "\n"
	     << "Lame lambda:       "   <<  lame_coefficient_lambda << "\n"
	     << std::endl;
   
   std::set<test_cases> test_cases_to_apply_pressurized_bc = 
     {test_cases::PRESSURIZED, test_cases::SNEDDON3D, test_cases::HET3D};
   std::set<test_cases> remaining_test_cases = current_test_case.set_diff(test_cases_to_apply_pressurized_bc);
   for (unsigned int i=0;i< pred_corr_levels;i++)
     {
       if (current_test_case.any(remaining_test_cases))
            {
              // TODO: constraints.close()?
              ConstraintMatrix constraints;
              constraints.close();
              
              std::vector<bool> component_mask (dim+1+dim, true);
              
              if (sub_test_case == "hohlwalzung")
                {
                  VectorTools::project (dof_handler,
                      constraints,
                      QGauss<dim>(degree+2),
                      InitialValuesHohlwalzung<dim>(alpha_eps),
                      solution
                      );
                }
                else if (current_test_case == test_cases::DYNAMIC_SLIT)
                  {
                  VectorTools::project (dof_handler,
                      constraints,
                      QGauss<dim>(degree+2),
                      InitialValuesDynamicSlit<dim>(),
                      solution
                      ); 
                  }
                else 
                {
                  VectorTools::project (dof_handler,
                      constraints,
                      QGauss<dim>(degree+2),
                      InitialValuesPhaseField<dim>(alpha_eps,
                                                    bool_initial_crack_via_phase_field,
                                                    current_test_case),
                      solution
                              );
                }
              
              compute_stress_per_cell();
              output_results (0,solution);
            }
       else if ( current_test_case.Is_pressurized() ||
		              current_test_case.Is_snedon())
      {
        /*
        ConstraintMatrix constraints;
        constraints.close();
        
        std::vector<bool> component_mask (dim+1+dim, true);
        VectorTools::project (dof_handler,
            constraints,
            QGauss<dim>(degree+2),
            InitialValuesPressurized<dim>(alpha_eps),
            solution
            );
        */
        
        VectorTools::interpolate(dof_handler,
                InitialValuesPressurized<dim>(alpha_eps), 
                solution);

        output_results (i,solution);
      }
       else if (current_test_case.Is_snedon3d())
          {
            ConstraintMatrix constraints;
            constraints.close();
            
            std::vector<bool> component_mask (dim+1+dim, true);
            VectorTools::project (dof_handler,
                constraints,
                QGauss<dim>(degree+2),
                InitialValuesSneddon3D<dim>(alpha_eps),
                solution
                );
            
            output_results (0,solution);
          }
       else if (current_test_case.Is_het3d())
          {
            ConstraintMatrix constraints;
            constraints.close();
            
            std::vector<bool> component_mask (dim+1+dim, true);
            VectorTools::project (dof_handler,
                constraints,
                QGauss<dim>(degree+2),
                InitialValuesHet3D<dim>(alpha_eps),
                solution
                );
            
            output_results (0,solution);
          }
       
       cout<<"PRE REFINEMENT STEP :: " << i << endl;
       refine_mesh();
      	 
       
     }

   // RESET THE INITIAL VALUES
   
   if (current_test_case.Is_pressurized() ||
       current_test_case.Is_snedon())
     {
        VectorTools::interpolate(dof_handler,
				    InitialValuesPressurized<dim>(alpha_eps), 
				    solution);

     }
   
   project_back_phase_field();

   
   
   if (bool_use_strain_history)
     {
       // Initialize strain history
       bool_set_initial_strain_history = true;
       update_quadrature_point_history ();
       output_results (0,solution);
       bool_set_initial_strain_history = false;
     }

  unsigned int refine_mesh_1st = 500000;
  unsigned int refine_mesh_2nd = 20000;
  unsigned int refine_mesh_3rd = 30000;
  unsigned int refine_mesh_4th = 40000;
  unsigned int refine_mesh_5th = 500000;
 
  const unsigned int output_skip = 1;
  unsigned int refinement_cycle = 0;

  double tmp_timestep;

  
  // Initialize old and old_old timestep sizes
  old_timestep = timestep;
  old_old_timestep = timestep;

  if (current_test_case.Is_pressurized())
    {
      current_pressure = 0.0; //1.0e-1; 
      old_timestep_current_pressure = 0.0; //1.0e-1; 
    }
  else 
    {
      current_pressure = 0.0;
      old_timestep_current_pressure = 0.0; //1.0e-1; 
    }

  double pressure_increment = 1.0e-1;

  double REL_functional_error = 0.0;
  double REL_old_timestep_functional_error = 0.0;

 // compute functional values special fields
  std::set<test_cases> speical_test_cases_functional_values = {test_cases::SNEDDON,test_cases::SNEDDON3D, test_cases::SCREW_DOMI};
  remaining_test_cases_functional_values = current_test_case.set_diff(speical_test_cases_functional_values);

  // Time loop
  do
    { 

      // Possibly adapt time step size
      //if (timestep_number > 5)
      //	timestep = 1.0e-2;

      //     if (timestep_number > 6)
      //	timestep = 1.0e-8;
//      if (timestep_number > 9)
//	{
//	  pressure_increment = 0.01;
//	  //gamma_penal = 1.0e+4;
//	}
//      else
//	{
//	  //pressure_increment = 1.0e+3;
//	  //gamma_penal = 1.0e+5;
//	}


      old_timestep_current_pressure = current_pressure;
      if (current_test_case.Is_pressurized())
	        current_pressure += pressure_increment; //5.0e-2 + time * 5.0e-2; //1.0e-1 + time * 1.0e-1;
      else if ((current_test_case.Is_snedon()) || ( current_test_case.Is_snedon3d()))
	        current_pressure = 1.0e-3;
      else if (current_test_case.Is_het3d())
	        current_pressure = 1.0e-3 + time * 0.25;
      else 
        	current_pressure = 0.0;

     

      if (current_test_case.Is_pressurized() || current_test_case.Is_het3d())
	      std::cout << "Current pressure: " << time << "\t" << current_pressure << std::endl;
      
      std::cout << std::endl;
      
      // Compute next time step
      total_number_newton_steps = 0;
      old_old_timestep_solution = old_timestep_solution;
      old_timestep_solution = solution;


      // Update time step sizes
      old_old_timestep = old_timestep;
      old_timestep = timestep;

      tmp_solution_for_time_adaptivity = solution;


      if (use_time_step_control)
          {
          if (timestep_number <= number_of_nonadaptive_time_steps) 
            time_stepping_scheme = "BE";
          else
            time_stepping_scheme = "BE_adaptive"; 

          std::cout << "Not tested in detail." << std::endl;
          abort();
          
        }

 
     if (time_stepping_scheme == "BE")
      {
        solve_spatial_problem ();
        time += timestep;
      }
     else if (time_stepping_scheme == "CN")
        {
    // n = 0 BE to obtain "smooth" solution
    if (timestep_number == 0)
      {
        theta = 1.0;
      }
    else
      theta = 0.5;
    
    solve_spatial_problem ();
    time += timestep;

    
        }
      else if (time_stepping_scheme == "BE_adaptive")
      {
	tmp_timestep = timestep;
	double tmp_time = time;
	unsigned int counter_redo_timestep = 0;

  

	// Adaptive time step control loop
	do {
	  if (counter_redo_timestep > 0)
	    {
	      std::cout << "Redo time step. Accuracy not good enough.\n" << std::endl;
	    }

	  solution = tmp_solution_for_time_adaptivity;
	  time = tmp_time;
	  tmp_timestep = timestep;

	  // Large time step 2k
	  timestep = 2.0 * timestep;
	  old_timestep_solution = tmp_solution_for_time_adaptivity;
	  solve_spatial_problem ();

	  // Evaluate goal functional
	  double gf_stress_x_2k = goal_functional_stress_x ();



	  std::cout << "-------------------------" << std::endl;
	  std::cout << "Two small time steps\n" << std::endl;
	  // Two times small time steps
	  time = tmp_time;
	  timestep = tmp_timestep;
	  solution = tmp_solution_for_time_adaptivity;

	  for (unsigned int l=0;l<2;l++)
	    {
	      old_timestep_solution = solution;
	      solve_spatial_problem ();

	       time += timestep; 
	    }

	  // Evaluate goal functional
	  double gf_stress_x_1k = goal_functional_stress_x ();


	  // Error evaluation
	  // ABS 
	  double ABS_functional_error = 0.0;

	  ABS_functional_error= std::abs(gf_stress_x_2k - gf_stress_x_1k);
	  REL_functional_error = ABS_functional_error / std::abs(gf_stress_x_1k);

	  double gamma_safety_factor = 1.0;
	  double power_time_step_control = 0.06666666666666666666;
	  double theta_time_step_control = gamma_safety_factor * 
	    std::pow((TOL_for_time_adaptivity) / REL_functional_error, power_time_step_control); 
	  
	  if (use_time_step_control)
	    {
	      if (theta_time_step_control < 1.0 ||
		  1.2 < theta_time_step_control)
		{
		  timestep = theta_time_step_control * timestep;
		}
	      else 
		timestep = 1.0 * timestep;
	    }
	  
	// PI-version (John and Rang) - does not yield satisfactory results
		//double rho_time_adapt = 0.9;
		//timestep = rho_time_adapt * timestep * timestep / old_timestep * 
		//  std::sqrt(TOL_for_time_adaptivity * REL_old_timestep_functional_error / (REL_functional_error * REL_functional_error));
		
		// Growth, maximum and minimum timesteps
		if (timestep > (timestep_growth_factor * tmp_timestep))
		  timestep = timestep_growth_factor * tmp_timestep;
		
		if (timestep >= timestep_upper_bound) 
		  timestep = timestep_upper_bound;
		
		if (timestep <= timestep_lower_bound)
		  timestep = timestep_lower_bound;
		
		std::cout << std::endl;
		std::cout << "-------------------------" << std::endl;
		std::cout << "Functional:   " << std::setprecision(5) << time << "   " <<  std::setprecision(10) << gf_stress_x_1k << std::endl;

		std::cout << "TOL:          " << std::setprecision(5) << time << "   " << std::setprecision(10) << TOL_for_time_adaptivity << std::endl;
		std::cout << "ABS error:    " << std::setprecision(5) << time << "   " << std::setprecision(10) << ABS_functional_error << std::endl;
		std::cout << "REL error:    " << std::setprecision(5) << time << "   " << std::setprecision(10) << REL_functional_error << std::endl;
		std::cout << "New timestep: " << std::setprecision(5) << time << "   " << std::setprecision(10) << timestep << std::endl;
		std::cout << "Old timestep: " << std::setprecision(5) << time << "   " << std::setprecision(10) << tmp_timestep << std::endl;
		std::cout << std::endl;
		
		counter_redo_timestep += 1;
		
	      } while (timestep < timestep_rejection_factor * tmp_timestep);
	      
	REL_old_timestep_functional_error = REL_functional_error;
	old_timestep = tmp_timestep;


      }


      // Inside adaptive time loop now
      //time += timestep;
      //timestep = tmp_timestep;
      
	
      // Compute functional values
      // Evaluate the summed goal functionals with the 'old' timestep,
      // thus we need to copy from tmp_timestep
      double new_timestep = timestep;


      // Using time step control, we always advance three time steps.
      // Thus we need to multiply with '2' for the goal functional evaluations.
      // TODO: introduce variable for time_adaptivity
      if (time_stepping_scheme == "BE_adaptive" && timestep_number > 0)
	timestep = 2.0 * tmp_timestep;
      else 
	timestep = 1.0 * tmp_timestep;



      // Compute functional values
      std::cout << std::endl;
      if (current_test_case.any(remaining_test_cases_functional_values))
      	compute_functional_values();

      if (current_test_case.Is_snedon() )
	      compute_functional_values_Sneddon();

      if (current_test_case.Is_screwdomi() )
        {
          compute_stress_per_cell();
          compute_energy();
        }

      if (current_test_case.Is_snedon3d() )
      	compute_cod_Sneddon3D();


      // Reset to standard timestep
      timestep = new_timestep;
      
      // Write solutions 
      // Plot corrector solution on new mesh (i.e., the final solution)
      if ((timestep_number % output_skip == 0))
	      output_results (timestep_number+1,solution);
      

      if (bool_use_strain_history)
	{
	  // Update stress history
	  update_quadrature_point_history ();
	}

      //      if (timestep_number  == refine_mesh_1st ||
      //  timestep_number  == refine_mesh_2nd ||
      //  timestep_number  == refine_mesh_3rd ||
      //  timestep_number  == refine_mesh_4th ||
      //  timestep_number  == refine_mesh_5th
      //  )			      			      			     
	{
	  std::cout << "Refinement cycle " 
		    << refinement_cycle 
		    << "\n================== "
		    << std::endl;
	  
	  refine_mesh ();
	  ++refinement_cycle;		
	}


 
      

      
      ++timestep_number;

    }
  while ((timestep_number <= max_no_timesteps) && (time <= end_time_value));
  
  
}

// The main function looks almost the same
// as in all other deal.II tuturial steps. 
int main () 
{
  try
    {
      deallog.depth_console (0);

      const unsigned int dimension = 2;
      Dynamic_Fracture_Problem<dimension> flow_problem(1);
      flow_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}




