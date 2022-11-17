#include "./../../include/dynamic_fracture.h"

// Here, we impose boundary conditions
// for the whole system.
template <int dim>
void
Dynamic_Fracture_Problem<dim>::set_initial_bc (const double time)
{ 
    std::map<unsigned int,double> boundary_values;  
    std::vector<bool> component_mask (dim+1+dim, false);
    component_mask[dim] = false; // scalar-valued phase-field
    component_mask[dim+1] = true;
    component_mask[dim+2] = true;
 
    if (test_case == "miehe_tension")
      {
    component_mask[0]     = false;
    component_mask[1]     = false;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      ZeroFunction<dim>(dim+1+dim),
					      //NonhomDirichletBoundaryValues<dim>(time),
					      boundary_values,
					      component_mask);  


    component_mask[0] = false;
    component_mask[1] = true;
    component_mask[2] = false; // phase_field
    VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
					      ZeroFunction<dim>(dim+1+dim),
					      //NonhomDirichletBoundaryValues2<dim>(time),
                                              boundary_values,
                                              component_mask);
 
    component_mask[0]   = false; // ux
    component_mask[1]   = true; 
    component_mask[2]   = false; // phase_field
    VectorTools::interpolate_boundary_values (dof_handler,
					      3,
					      NonhomDirichletBoundaryValues<dim>(time, test_case,alpha_eps),
					      boundary_values,
					      component_mask);

		//    component_mask[0] = false;
		//    component_mask[1] = false;
		//    component_mask[2] = true; // phase_field
		//    VectorTools::interpolate_boundary_values (dof_handler,
		//                                              2,
		//					      ConstantFunction<dim>(1.0, dim+1),
		//                                              boundary_values,
		//                                              component_mask);
		//
		//    component_mask[0]   = false; // ux
		//    component_mask[1]   = false;
		//    component_mask[2]   = true; // phase_field
		//    VectorTools::interpolate_boundary_values (dof_handler,
		//					      3,
		//					      ConstantFunction<dim>(1.0, dim+1),
		//					      boundary_values,
		//					      component_mask);
	}
	else if (test_case == "miehe_shear")
	{
		// ux = component_mask[0]
		// uy = component_mask[1]
		// phi = component_mask[2]
		// vx = component_mask[3]
		// vy = component_mask[4]
		component_mask[0] = false;
		component_mask[1] = true;

		// changed component_mask[dim+1]     = false;
		component_mask[dim + 1] = false;
		component_mask[dim + 2] = true;
		// 0 is the left edge
		VectorTools::interpolate_boundary_values(dof_handler,
												 0,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		// 1 is the right edge
		component_mask[0] = false;
		component_mask[1] = true;
		// changed component_mask[dim+1]     = false;
		component_mask[dim + 1] = false;
		component_mask[dim + 2] = true;

		VectorTools::interpolate_boundary_values(dof_handler,
												 1,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		// 2 is the bottom edge
		component_mask[0] = true;
		component_mask[1] = true;
		component_mask[dim + 1] = true;
		component_mask[dim + 2] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 2,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		// 3 is the top edge
		component_mask[0] = true;
		// changed component_mask[1]     = true;
		component_mask[1] = true;
		component_mask[dim + 1] = true;
		// changed component_mask[dim+2]     = true;
		component_mask[dim + 2] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 3,
												 NonhomDirichletBoundaryValues<dim>(time, test_case, alpha_eps),
												 boundary_values,
												 component_mask);

		// 4 is the bottom part of the crack
		// uy=0, vy=0
		component_mask[0] = false;
		component_mask[1] = true;
		component_mask[dim + 1] = false;
		component_mask[dim + 2] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 4,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		if (bool_initial_crack_via_phase_field)
		{
			// Phase-field
			component_mask[0] = false;
			component_mask[1] = false;
			component_mask[2] = true;
			component_mask[dim + 1] = false;
			component_mask[dim + 2] = false;
			VectorTools::interpolate_boundary_values(dof_handler,
													 1,
													 NonhomDirichletBoundaryValues<dim>(time, test_case, alpha_eps),
													 boundary_values,
													 component_mask);
		}
	}
	else if (test_case == "l_shaped")
	{
		component_mask[0] = false;
		component_mask[1] = false;
		VectorTools::interpolate_boundary_values(dof_handler,
												 0,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = false;
		component_mask[1] = false;
		VectorTools::interpolate_boundary_values(dof_handler,
												 11,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = true;
		component_mask[1] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 2,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = false;
		component_mask[1] = false;
		VectorTools::interpolate_boundary_values(dof_handler,
												 3,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = false;
		component_mask[1] = false;
		VectorTools::interpolate_boundary_values(dof_handler,
												 4,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = true;
		component_mask[1] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 5,
												 NonhomDirichletBoundaryValues<dim>(time, test_case, alpha_eps),
												 boundary_values,
												 component_mask);

		component_mask[0] = false;
		component_mask[1] = false;
		component_mask[2] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 1,
												 // ZeroFunction<dim>(dim+1),
												 ConstantFunction<dim>(1.0, dim + 1 + dim),
												 boundary_values,
												 component_mask);

		VectorTools::interpolate_boundary_values(dof_handler,
												 5,
												 ConstantFunction<dim>(1.0, dim + 1 + dim),
												 boundary_values,
												 component_mask);
	}
	else if (test_case == "pressurized" || test_case == "Sneddon")
	{
		component_mask[0] = true;
		component_mask[1] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 0,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);
	}
	else if (test_case == "screw_domi")
	{
		component_mask[0] = false;
		component_mask[1] = false;
		VectorTools::interpolate_boundary_values(dof_handler,
												 0,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = false;
		component_mask[1] = false;
		VectorTools::interpolate_boundary_values(dof_handler,
												 1,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = false;
		component_mask[1] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 2,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		// TODO
		component_mask[0] = false; // false
		component_mask[1] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 3,
												 NonhomDirichletBoundaryValues<dim>(time, test_case, alpha_eps),
												 boundary_values,
												 component_mask);

		component_mask[0] = false;
		component_mask[1] = false;
		component_mask[2] = false;
		VectorTools::interpolate_boundary_values(dof_handler,
												 2,
												 // ZeroFunction<dim>(dim+1),
												 ConstantFunction<dim>(1.0, dim + 1 + dim),
												 boundary_values,
												 component_mask);
	}
	else if ((test_case == "Sneddon3D") || (test_case == "Het3D"))
	{
		component_mask[0] = true;
		component_mask[1] = true;
		component_mask[2] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 0,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = true;
		component_mask[1] = true;
		component_mask[2] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 1,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = true;
		component_mask[1] = true;
		component_mask[2] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 2,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = true;
		component_mask[1] = true;
		component_mask[2] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 3,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = true;
		component_mask[1] = true;
		component_mask[2] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 4,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);

		component_mask[0] = true;
		component_mask[1] = true;
		component_mask[2] = true;
		VectorTools::interpolate_boundary_values(dof_handler,
												 5,
												 ZeroFunction<dim>(dim + 1 + dim),
												 boundary_values,
												 component_mask);
	}

	for (typename std::map<unsigned int, double>::const_iterator
			 i = boundary_values.begin();
		 i != boundary_values.end();
		 ++i)
		solution(i->first) = i->second;

}
