#include "./../../include/dynamic_fracture.h"

template<int dim>
bool Dynamic_Fracture_Problem<dim>::refine_mesh()
{

  std::string refinement_type = "pff_refinement";
  //std::string refinement_type= "fixed_refinement_2D";
  //std::string refinement_type = "fixed_refinement_Sneddon3D";


  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();


  if (refinement_type == "fixed_refinement_2D")
    {
      for (; cell!=endc; ++cell)
	{
	  
	  for (unsigned int vertex=0;vertex < GeometryInfo<dim>::vertices_per_cell;++vertex)	   	 
	    {
	      Tensor<1,dim> cell_vertex = (cell->vertex(vertex));
	      if (cell_vertex[0] <= 0.55 && cell_vertex[0] >= 0.0 && 
		  cell_vertex[1] < 0.51 && cell_vertex[1] > -0.1) 	  
		{
		  cell->set_refine_flag();
		  break;
		}
	    }
	}

    }
  else if (current_test_case.Is_snedon3d())
    {
      for (; cell!=endc; ++cell)
	{
	  	  std::cout << "Drin" << std::endl;
	  for (unsigned int vertex=0;vertex < GeometryInfo<dim>::vertices_per_cell;++vertex)	   	 
	    {
	      Tensor<1,dim> cell_vertex = (cell->vertex(vertex));
	      if (cell_vertex[0] <= 6.5 && cell_vertex[0] >= 3.5 && 
		  cell_vertex[2] <= 6.5 && cell_vertex[2] >= 3.5 && 
		  cell_vertex[1] <= 5.5 && cell_vertex[1] >= 4.5) 	  
		{
		  cell->set_refine_flag();
		  break;
		}
	    }
	  
	}
      

    }
  else if (current_test_case.Is_het3d())
    {
      for (; cell!=endc; ++cell)
	{

	  for (unsigned int vertex=0;vertex < GeometryInfo<dim>::vertices_per_cell;++vertex)	   	 
	    {
	      Tensor<1,dim> cell_vertex = (cell->vertex(vertex));
	      if (cell_vertex[0] <= 7.0 && cell_vertex[0] >= 2.0 && 
		  cell_vertex[2] <= 7.0 && cell_vertex[2] >= 3.0 && 
		  cell_vertex[1] <= 7.0 && cell_vertex[1] >= 3.0) 	  
		{
		  cell->set_refine_flag();
		  break;
		}
	    }
	  
	}
      

    }
  else if (refinement_type == "pff_refinement")
    {
      std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);

 for (; cell!=endc; ++cell)
    {

     cell->get_dof_indices(local_dof_indices);
     for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	{
	  const unsigned int comp_i = fe.system_to_component_index(i).first;
	  if (comp_i != dim)
	    continue; // only look at phase field
	  if (solution(local_dof_indices[i])
	      < value_phase_field_for_refinement )
	    {
	      cell->set_refine_flag();
	      break;
	    }
	}

      // Global refinement
      //cell->set_refine_flag();

    }


 unsigned int total_refinement_levels = global_refinement_steps + pred_corr_levels;
 // Limit number of refinement levels
 // if (test_case != TestCase::sneddon_2d)
    {
      typename DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();
 
     for (; cell != endc; ++cell)
       {
	 if (cell->level() == static_cast<int>(total_refinement_levels))
	   cell->clear_refine_flag();
       }
    }

    } // end pff_refinement



 // Check if cells are flagged for refinement 
 // or coarsening
 {
   bool refine_or_coarsen = false;
   triangulation.prepare_coarsening_and_refinement();
    
    typename DoFHandler<dim>::active_cell_iterator cell =
      dof_handler.begin_active(), endc = dof_handler.end();

    for (; cell != endc; ++cell)
      if (cell->refine_flag_set() || cell->coarsen_flag_set())
        {
          refine_or_coarsen = true;
          break;
        }

    if (refine_or_coarsen == 0)
      return false;
  }


  BlockVector<double> tmp_solution;
  tmp_solution = solution;

  BlockVector<double> tmp_old_timestep_solution;
  tmp_old_timestep_solution = old_timestep_solution;

  BlockVector<double> local_tmp_solution;
  local_tmp_solution = tmp_solution_for_time_adaptivity;

  SolutionTransfer<dim, BlockVector<double> > solution_transfer (dof_handler);
  SolutionTransfer<dim, BlockVector<double> > solution_transfer_2 (dof_handler);
  SolutionTransfer<dim, BlockVector<double> > solution_transfer_3 (dof_handler);
  
  //triangulation.prepare_coarsening_and_refinement();
  solution_transfer.prepare_for_coarsening_and_refinement(tmp_solution);
  solution_transfer_2.prepare_for_coarsening_and_refinement(tmp_old_timestep_solution);
  solution_transfer_3.prepare_for_coarsening_and_refinement(local_tmp_solution);
  
  triangulation.execute_coarsening_and_refinement ();
  setup_system ();

  solution_transfer.interpolate(tmp_solution, solution);
  solution_transfer_2.interpolate(tmp_old_timestep_solution, old_timestep_solution); 
  solution_transfer_3.interpolate(local_tmp_solution, tmp_solution_for_time_adaptivity); 

  return true;

}

