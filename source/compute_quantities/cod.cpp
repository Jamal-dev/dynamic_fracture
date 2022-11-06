#include "./../../include/dynamic_fracture.h"

template <int dim>
double
Dynamic_Fracture_Problem<dim>::compute_cod (
  const double eval_line)
{

  const QGauss<dim - 1> face_quadrature_formula(3);
  FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
                                   update_values | update_quadrature_points | update_gradients
                                   | update_normal_vectors | update_JxW_values);



  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_face_q_points = face_quadrature_formula.size();

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  std::vector<Vector<double> > face_solution_values(n_face_q_points,
                                                    Vector<double>(dim+1+dim));
  std::vector<std::vector<Tensor<1, dim> > > face_solution_grads(
    n_face_q_points, std::vector<Tensor<1, dim> >(dim+1+dim));

  double cod_value = 0.0;
  double eps = 1.0e-6;

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active(), endc = dof_handler.end();

  for (; cell != endc; ++cell)
      {
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          {
            fe_face_values.reinit(cell, face);
            fe_face_values.get_function_values(solution,
                                               face_solution_values);
            fe_face_values.get_function_gradients(solution,
                                                  face_solution_grads);

            for (unsigned int q_point = 0; q_point < n_face_q_points;
                 ++q_point)
              {
                if ((fe_face_values.quadrature_point(q_point)[0]
                     < (eval_line + eps))
                    && (fe_face_values.quadrature_point(q_point)[0]
                        > (eval_line - eps)))
                  {
                    const Tensor<1, dim> u = ALE_Transformations::get_u<dim>(
                                               q_point, face_solution_values);

                    const Tensor<1, dim> grad_pf =
                      ALE_Transformations::get_grad_pf<dim>(q_point,
                                                face_solution_grads);

                    cod_value += 0.25 * u * grad_pf
                                 * fe_face_values.JxW(q_point);

                  }

              }
          }
      }


  std::cout << eval_line << "  " << cod_value << std::endl;

  return cod_value;

}

