#include "./../../include/dynamic_fracture.h"


// Now, we arrive at the function that is responsible 
// to compute the line integrals for stress evaluations
template <int dim>
double Dynamic_Fracture_Problem<dim>::goal_functional_stress_x ()
{
    
  const QGauss<dim-1> face_quadrature_formula (3);
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
				    update_values | update_gradients | update_normal_vectors | 
				    update_JxW_values);
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_face_q_points = face_quadrature_formula.size();

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  std::vector<Vector<double> >  face_solution_values (n_face_q_points, 
						      Vector<double> (dim+1+dim));

  std::vector<std::vector<Tensor<1,dim> > > 
    face_solution_grads (n_face_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));
  
  Tensor<1,dim> load_value;

 const Tensor<2, dim> Identity =
    ALE_Transformations::get_Identity<dim>();

 unsigned int bc_color = 3;
 if (current_test_case.Is_miehe_shear())
   bc_color = 3;
 else if (current_test_case.Is_l_shaped())
   bc_color = 5;
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  unsigned int cell_counter = 0;
   for (; cell!=endc; ++cell)
     {
       //lame_coefficient_mu = lame_coefficient_mu_vector(cell_counter);
       //lame_coefficient_lambda = lame_coefficient_lambda_vector(cell_counter);
       cell_counter++;



       for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	 if (cell->face(face)->at_boundary() && 
	     // TODO: 3 or 5
	     cell->face(face)->boundary_id() == bc_color)
	   { 
	     fe_face_values.reinit (cell, face);
	     fe_face_values.get_function_values (solution, face_solution_values);
	     fe_face_values.get_function_gradients (solution, face_solution_grads);
	 	      
	     for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
	       {
		 const Tensor<2,dim> grad_u = ALE_Transformations 
		   ::get_grad_u<dim> (q_point, face_solution_grads);
		 
		 const Tensor<2, dim> E = 0.5 * (grad_u + transpose(grad_u));
		 const double tr_E = dealii::trace(E);
		 
		 Tensor<2, dim> stress_term;
		 stress_term = lame_coefficient_lambda * tr_E * Identity
		   + 2 * lame_coefficient_mu * E;
		 
		 load_value +=  stress_term *
		   fe_face_values.normal_vector(q_point)* fe_face_values.JxW(q_point);

		 
	       }
	   } // end boundary 3 for structure
       

     } 

    // load_cases_started
   if (current_test_case.Is_miehe_shear()  || 
			 current_test_case.Is_p_notched_cavity() ||
       current_test_case.Is_miehe_tension()
       || current_test_case.Is_pmesh1()) // load_cases_ended
     {
       load_value[0] *= -1.0;


     }
   else if (current_test_case.Is_l_shaped())
     {
       // TODO (change depending on the mesh)
       load_value[1] *= -1.0;

       double load_increment = time;
        if (time < 0.3)
          load_increment = time;
        else if (time >= 0.3 && time < 0.8)
          load_increment = 0.6 - time;
        else if (time >= 0.8)
          load_increment = -1.0 + time;


     }

   return load_value[0];

}

