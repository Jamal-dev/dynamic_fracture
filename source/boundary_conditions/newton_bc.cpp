#include "./../../include/dynamic_fracture.h"

// This function applies boundary conditions 
// to the Newton iteration steps. For all variables that
// have Dirichlet conditions on some (or all) parts
// of the outer boundary, we apply zero-Dirichlet
// conditions, now. 
template <int dim>
void
Dynamic_Fracture_Problem<dim>::set_newton_bc ()
{
    // ConstraintMatrix constraints; // let it use the global constraint
	std::map<unsigned int,double> boundary_values; 
	std::vector<bool> component_mask (dim+1+dim, false);
    component_mask[dim] = false;  // false  // phase_field
    component_mask[dim+1] = true;
    component_mask[dim+2] = true;

    if (current_test_case == test_cases::MIEHE_TENSION)
      {
			component_mask[0]       = false;
			component_mask[1]     = false;
			component_mask[dim+1]     = false;
			component_mask[dim+2]     = false;
			VectorTools::interpolate_boundary_values (dof_handler,
								0,
								ZeroFunction<dim>(dim+1+dim), 
								constraints,				
								component_mask);  

		
			component_mask[0] = false;
			component_mask[1] = true;
			component_mask[2] = false;
			component_mask[dim+1]     = false;
			component_mask[dim+2]     = false;
			VectorTools::interpolate_boundary_values (dof_handler,
													2,
								ZeroFunction<dim>(dim+1+dim),  
								constraints,
													component_mask);
		
			component_mask[0] = false;
			component_mask[1] = true; 
			component_mask[2] = false;
			component_mask[dim+1]     = false;
			component_mask[dim+2]     = false;
			VectorTools::interpolate_boundary_values (dof_handler,
								3,
								ZeroFunction<dim>(dim+1+dim),  
								constraints,		
								component_mask);

      }
	else if (current_test_case == test_cases::P_MESH_1)
      {
		/*
			left_edge,      1
			right_edge,     0
			top_edge,       3
			bottom_edge,    2
			crack_bottom,   4
			crack_top,      5
		*/
		component_mask[0]       = false;
		component_mask[1]     = false;
		component_mask[dim+1]     = false;
		component_mask[dim+2]     = false;
		VectorTools::interpolate_boundary_values (dof_handler,
							0,
							ZeroFunction<dim>(dim+1+dim), 
							constraints,				
							component_mask);  

	
		component_mask[0] = false;
		component_mask[1] = true;
		component_mask[2] = false;
		component_mask[dim+1]     = false;
		component_mask[dim+2]     = false;
		VectorTools::interpolate_boundary_values (dof_handler,
												2,
							ZeroFunction<dim>(dim+1+dim),  
							constraints,
												component_mask);
	
		// adding condition on both uy and vy
		component_mask[0] = false;
		component_mask[1] = true; 
		component_mask[2] = false;
		component_mask[dim+1]     = false;
		component_mask[dim+2]     = true;
		VectorTools::interpolate_boundary_values (dof_handler,
							3,
							ZeroFunction<dim>(dim+1+dim),  
							constraints,		
							component_mask);

	}
else if (current_test_case == test_cases::P_ASYMMETRY)
	{
		/*
			
		*/
			component_mask[0]     = false;
			component_mask[1]     = true;
			component_mask[dim+1]     = false;
			component_mask[dim+2]     = false;
			VectorTools::interpolate_boundary_values (dof_handler,
								3,
								ZeroFunction<dim>(dim+1+dim),
								boundary_values,
								component_mask);
			component_mask[0]     = false;
			component_mask[1]     = false;
			component_mask[dim+1]     = false;
			component_mask[dim+2]     = true;
			VectorTools::interpolate_boundary_values (dof_handler,
								3,
								ZeroFunction<dim>(dim+1+dim),
								boundary_values,
								component_mask);  

			/*
			  Left and right corners fixing
			*/
			// fix y component of left and right bottom corners
			typename DoFHandler<dim>::active_cell_iterator cell =
				dof_handler.begin_active(), endc = dof_handler.end();

			for (; cell != endc; ++cell)
				{
				if (cell->is_artificial())
					continue;

				for (unsigned int v = 0;
					v < GeometryInfo<dim>::vertices_per_cell; ++v)
					{
					if (
						std::abs(cell->vertex(v)[1]) < 1e-10
						&&
						(
						std::abs(cell->vertex(v)[0]-1.0) < 1e-10
						|| std::abs(cell->vertex(v)[0]-19.0) < 1e-10
						))
						{
						types::global_dof_index idx = cell->vertex_dof_index(v, 1);// 1=y displacement
						boundary_values[idx] = 0.0;


						idx = cell->vertex_dof_index(v, 0);// 0=x displacement
						if (std::abs(cell->vertex(v)[0]-1.0) < 1e-10)
							boundary_values[idx] = 0.0;
				// Info: We only fix the x-value of the left bc
				//else if (std::abs(cell->vertex(v)[0]-19.0) < 1e-10)
				//  boundary_values[idx] = 0.0;

						//idx = cell->vertex_dof_index(v, 2);// 2= phase-field
						//boundary_values[idx] = 1.0;
						}
				// pushing on top
		//              else if (
		//                std::abs(cell->vertex(v)[0]-10.0) < 1e-10
		//                &&
		//                std::abs(cell->vertex(v)[1]-8.0) < 1e-10
		//              )
		//                {
		//                  types::global_dof_index idx = cell->vertex_dof_index(v, 1);// 1=y displacement
		//                  boundary_values[idx] = -1.0*time;
		//                }
		//
					} // end running of vertices per cell
				} // end cells

			/*
			  corner fixing finished
			*/
		
		//    component_mask[0]     = true; 
		//    component_mask[1]     = true; 
		//    component_mask[3]     = true; 
		//    component_mask[4]     = true; 
		//    VectorTools::interpolate_boundary_values (dof_handler,
		//					      40,
		//					      ZeroFunction<dim>(dim+1+dim),  
		//					      constraints,		
		//					      component_mask);
		//
		//    component_mask[0]     = true; 
		//    component_mask[1]     = true;
		//    component_mask[3]     = true; 
		//    component_mask[4]     = true; 
		//    VectorTools::interpolate_boundary_values (dof_handler,
		//					      41,
		//					      ZeroFunction<dim>(dim+1+dim),  
		//					      constraints,		
		//					      component_mask);
		//
		//    component_mask[0]     = true; 
		//    component_mask[1]     = true;
		//    component_mask[3]     = true; 
		//    component_mask[4]     = true; 
		//    VectorTools::interpolate_boundary_values (dof_handler,
		//					      42,
		//					      ZeroFunction<dim>(dim+1+dim),  
		//					      constraints,		
		//					      component_mask);
		//
			


	} //end case asymmetry
else if (current_test_case == test_cases::P_NOTCHED_CAVITY)
	{
		/*
			left_edge,      1
			right_edge,     0
			top_edge,       3
			bottom_edge,    2
			circle_edge,    5
		*/
		// fix hole to 0
		component_mask[0]     = true;
		component_mask[1]     = true;
		component_mask[2]     = false;
		component_mask[dim+1]     = true;
		component_mask[dim+2]     = true;
		VectorTools::interpolate_boundary_values (dof_handler,
							5,
							ZeroFunction<dim>(dim+1+dim),
							//NonhomDirichletBoundaryValues<dim>(time),
							boundary_values,
							component_mask);  

	
		component_mask[0] = false;
		component_mask[1] = true;
		component_mask[2] = false;
		component_mask[dim+1]     = false;
		component_mask[dim+2]     = false;
		VectorTools::interpolate_boundary_values (dof_handler,
												2,
							ZeroFunction<dim>(dim+1+dim),  
							constraints,
												component_mask);
	
		// adding condition on both uy and vy
		component_mask[0] = false;
		component_mask[1] = true; 
		component_mask[2] = false;
		component_mask[dim+1]     = false;
		component_mask[dim+2]     = true;
		VectorTools::interpolate_boundary_values (dof_handler,
							3,
							ZeroFunction<dim>(dim+1+dim),  
							constraints,		
							component_mask);

}

	else if  (current_test_case == test_cases::DYNAMIC_SLIT)
      {
		// ux = component_mask[0]
        // uy = component_mask[1]
        // phi = component_mask[2]
        // vx = component_mask[3]
        // vy = component_mask[4]

		// 0 is the left edge
		// 1 is the right edge
		// 2 is the bottom edge
		// 3 is the top edge
		// 4 is the crack line


    component_mask[0] = false;
    component_mask[1] = true;
    component_mask[2] = false; // phase_field
	component_mask[dim+1]     = false;
    component_mask[dim+dim]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              2/*bottom edge*/,
					      ZeroFunction<dim>(dim+1+dim),
					      //NonhomDirichletBoundaryValues2<dim>(time),
                                              constraints,
                                              component_mask);
 
    component_mask[0]   = false; // ux
    component_mask[1]   = true; 
    component_mask[2]   = false; // phase_field
	component_mask[dim+1]     = false;
    component_mask[dim+dim]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
					      3/*top edge*/,
					      ZeroFunction<dim>(dim+1+dim),
					      constraints,
					      component_mask);


      }
    else if (current_test_case == test_cases::MIEHE_SHEAR)
      {
		/*
			left_edge,      0
			right_edge,     1
			top_edge,       3
			bottom_edge,    2
			crack_bottom,   4
			crack_top,      5
		*/
        // 0 is the left edge
	component_mask[0]     = false;
	component_mask[1]     = true;
	component_mask[dim+1]     = false;
	component_mask[dim+2]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  0,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);  

	
  // 1 is the right edge
  component_mask[0]       = false;
	component_mask[1]     = true;
	component_mask[dim+1]     = false;
	component_mask[dim+2]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  1,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);

	// 2 is the bottom edge
  component_mask[0]     = true;
	component_mask[1]     = true;
	component_mask[dim+1]     = true;
	component_mask[dim+2]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  2,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);  

	
  // 3 is the top edge
  component_mask[0]     = true;
  // changed component_mask[1]     = true;
	component_mask[1]     = true;
	component_mask[dim+1]     = true;
	// changed component_mask[dim+2]     = true;
  component_mask[dim+2]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  3,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);  

	//4 is the crack line
  component_mask[0]     = false;
	component_mask[1]     = true;
	component_mask[dim+1]     = false;
	component_mask[dim+2]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  4,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);  

	if (bool_initial_crack_via_phase_field)
	{
	component_mask[0]     = false;
	component_mask[1]     = false;
	component_mask[2]     = true;
	component_mask[dim+1]     = false;
	component_mask[dim+2]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  1,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);
}


      }
    else if (current_test_case == test_cases::L_SHAPED)
      {
	component_mask[0]     = false;
	component_mask[1]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  0,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask); 



	component_mask[0]     = false;
	component_mask[1]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  11,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);

	component_mask[0]     = true;
	component_mask[1]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  2,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);  

	component_mask[0]     = false;
	component_mask[1]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  3,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);  

	component_mask[0]     = false;
	component_mask[1]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  4,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);

	component_mask[0]     = true;
	component_mask[1]     = true;
	component_mask[2]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  5,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);

	component_mask[0]     = false;
	component_mask[1]     = false;
	component_mask[2]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  1,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask); 

  

      }
    else if (current_test_case == test_cases::PRESSURIZED || current_test_case == test_cases::SNEDDON)
      {
	component_mask[0]     = true;
	component_mask[1]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  0,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);  

      }
    else if (current_test_case == test_cases::SCREW_DOMI)
      {
	component_mask[0]       = false;
	component_mask[1]       = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  0,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,		
						  component_mask);  
	
	
	component_mask[0]     = false;
	component_mask[1]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  1,
						  ZeroFunction<dim>(dim+1+dim),  
						  constraints,
						  component_mask);
	component_mask[0]     = false; // false
	component_mask[1]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  2,
						  ZeroFunction<dim>(dim+1+dim),  
						  constraints,
						  component_mask);
	
	// TODO - top boundary
	component_mask[0]     = false; // false
	component_mask[1]     = true; 
	VectorTools::interpolate_boundary_values (dof_handler,
						  3,
						  ZeroFunction<dim>(dim+1+dim),  
						  constraints,		
						  component_mask);

	// Phase-field bc values
	component_mask[0]     = false;
	component_mask[1]     = false;
	component_mask[2]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  2,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask); 

	

      }
    else if ((current_test_case == test_cases::SNEDDON3D) || (current_test_case == test_cases::HET3D))
      {
    component_mask[0]       = true;
    component_mask[1]       = true;
    component_mask[2]       = true;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      ZeroFunction<dim>(dim+1+dim), 
					      constraints,		
					      component_mask);  

 
    component_mask[0]     = true;
    component_mask[1]     = true;
    component_mask[2]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
					      1,
					      ZeroFunction<dim>(dim+1+dim),  
					      constraints,
					      component_mask);
    component_mask[0]     = true; // false
    component_mask[1]     = true;
    component_mask[2]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
					      ZeroFunction<dim>(dim+1+dim),  
					      constraints,
                                              component_mask);
 
    component_mask[0]     = true; // false
    component_mask[1]     = true;
    component_mask[2]     = true; 
    VectorTools::interpolate_boundary_values (dof_handler,
					      3,
					      ZeroFunction<dim>(dim+1+dim),  
					      constraints,		
					      component_mask);

  component_mask[0]     = true; // false
    component_mask[1]     = true;
    component_mask[2]     = true; 
    VectorTools::interpolate_boundary_values (dof_handler,
					      4,
					      ZeroFunction<dim>(dim+1+dim),  
					      constraints,		
					      component_mask);

    component_mask[0]     = true; // false
    component_mask[1]     = true;
    component_mask[2]     = true; 
    VectorTools::interpolate_boundary_values (dof_handler,
					      5,
					      ZeroFunction<dim>(dim+1+dim),  
					      constraints,		
					      component_mask);
      }


}  

