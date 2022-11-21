#include "./../../include/dynamic_fracture.h"

// This function is known from almost all other 
// tutorial steps in deal.II.
template <int dim>
void
Dynamic_Fracture_Problem<dim>::output_results (const unsigned int refinement_cycle,
			      const BlockVector<double> output_vector)  const
{

  std::vector<std::string> solution_names;
  solution_names.push_back ("x_dis");
  solution_names.push_back ("y_dis");
  if (dim == 3)
    solution_names.push_back ("z_dis");
  solution_names.push_back ("phase_field");

  solution_names.push_back ("x_vel");
  solution_names.push_back ("y_vel");
  if (dim == 3)
    solution_names.push_back ("z_vel");

 std::vector<std::string> solution_names_lambda_penal;
  solution_names_lambda_penal.push_back ("x_lp");
  solution_names_lambda_penal.push_back ("y_lp");
  if (dim == 3)
    solution_names_lambda_penal.push_back ("z_lp");

  solution_names_lambda_penal.push_back ("pf_lp");

  solution_names_lambda_penal.push_back ("vx_lp");
  solution_names_lambda_penal.push_back ("vy_lp");
  if (dim == 3)
    solution_names_lambda_penal.push_back ("vz_lp");
      

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim+1+dim, DataComponentInterpretation::component_is_scalar);


  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);  
   
  data_out.add_data_vector (output_vector, solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  
  // only required for heterogeneous materials
  //data_out.add_data_vector(lame_coefficient_mu_vector, "mu");
  //data_out.add_data_vector(lame_coefficient_lambda_vector, "lambda");

  if (current_test_case.IsScrewDomi())
    data_out.add_data_vector(solution_stress_per_cell, "stress"); 


  data_out.add_data_vector (solution_lambda_penal_func, solution_names_lambda_penal,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);


  Vector<float> e_mod(triangulation.n_active_cells());
 if (current_test_case.IsHet3D())
    {
 
      typename DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();

      unsigned int cellindex = 0;
      for (; cell != endc; ++cell, ++cellindex)
          {
            e_mod(cellindex) = 1.0 + func_emodulus->value(cell->center(), 0);
          }
      data_out.add_data_vector(e_mod, "E_modulus");
    }


  data_out.build_patches ();

  // Output
  // 1X XXXX := predictor output
  // 2X XXXX := corrector output
  // 3X XXXX := old solution on new mesh
  
  // X0 XXXX := level 0 of pred-corr
  // X1 XXXX := level 1 of pred-corr
  // etc.
   
  std::ostringstream filename;

  std::cout << "------------------" << std::endl;
  std::cout << "Write solution" << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << std::endl;
  filename << filename_basis
	   << Utilities::int_to_string (refinement_cycle, 6)
	   << ".vtk";
  
  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);

}


