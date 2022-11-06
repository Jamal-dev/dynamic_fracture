#include "./../../include/dynamic_fracture.h"

  template <int dim>
  void Dynamic_Fracture_Problem<dim>::setup_quadrature_point_history ()
  {
    QGauss<dim>   quadrature_formula(degree+2);

    unsigned int our_cells = 0;
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
        ++our_cells;

    triangulation.clear_user_data();


    {
      std::vector<PointHistory<dim> > tmp;
      tmp.swap (quadrature_point_history);
    }
    quadrature_point_history.resize (our_cells *
                                     quadrature_formula.size());


    unsigned int history_index = 0;
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      {
	cell->set_user_pointer (&quadrature_point_history[history_index]);
	history_index += quadrature_formula.size();
      }


    Assert (history_index == quadrature_point_history.size(),
            ExcInternalError());
  }

