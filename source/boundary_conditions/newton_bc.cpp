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
	std::vector<bool> component_mask (dim+1+dim, false);
    component_mask[dim] = false;  // false  // phase_field
    component_mask[dim+1] = true;
    component_mask[dim+2] = true;

    if (test_case == "miehe_tension")
      {
    component_mask[0]       = false;
    component_mask[1]     = false;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      ZeroFunction<dim>(dim+1+dim), 
					      constraints,				
					      component_mask);  

 
    component_mask[0] = false;
    component_mask[1] = true;
    component_mask[2] = false;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
					      ZeroFunction<dim>(dim+1+dim),  
					      constraints,
                                              component_mask);
 
    component_mask[0] = false;
    component_mask[1] = true; 
    component_mask[2] = false;
    VectorTools::interpolate_boundary_values (dof_handler,
					      3,
					      ZeroFunction<dim>(dim+1+dim),  
					      constraints,		
					      component_mask);

      }
	else if  (test_case == "dynamic_slit")
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
    else if (test_case == "miehe_shear")
      {
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
    else if (test_case == "l_shaped")
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
    else if (test_case == "pressurized" || test_case == "Sneddon")
      {
	component_mask[0]     = true;
	component_mask[1]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  0,
						  ZeroFunction<dim>(dim+1+dim), 
						  constraints,				
						  component_mask);  

      }
    else if (test_case == "screw_domi")
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
    else if ((test_case == "Sneddon3D") || (test_case == "Het3D"))
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

