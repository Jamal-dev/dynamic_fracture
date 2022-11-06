#include "./../../include/dynamic_fracture.h"

template <int dim>
void Dynamic_Fracture_Problem<dim>::compute_cod_Sneddon3D ()
{
  double x1,y1, y2, x2, x3, x4, x5, x6, x7, x8, x9,x10;  
  double y_and_eps = 5.0 + alpha_eps/2.0;
  x1 = 2.0 * compute_point_value(Point<dim>(4.0,   y_and_eps, 5.0), 1);
  x2 = 2.0 * compute_point_value(Point<dim>(4.25,  y_and_eps, 5.0), 1);
  x3 = 2.0 * compute_point_value(Point<dim>(4.5,   y_and_eps, 5.0), 1);
  x4 = 2.0 * compute_point_value(Point<dim>(4.75,  y_and_eps, 5.0), 1);
  x5 = 2.0 * compute_point_value(Point<dim>(5.0,   y_and_eps, 5.0), 1);
  x6 = 2.0 * compute_point_value(Point<dim>(5.25,  y_and_eps, 5.0), 1);
  x7 = 2.0 * compute_point_value(Point<dim>(5.5,   y_and_eps, 5.0), 1);
  x8 = 2.0 * compute_point_value(Point<dim>(5.75,  y_and_eps, 5.0), 1);
  x9 = 2.0 * compute_point_value(Point<dim>(6.0,   y_and_eps, 5.0), 1);

  std::cout << "------------------" << std::endl;
  std::cout << "DisY1:   " << "4.0  " << x1 << std::endl;
  std::cout << "DisY2:   " << "4.25 " << x2 << std::endl;
  std::cout << "DisY3:   " << "4.5  " << x3 << std::endl;
  std::cout << "DisY4:   " << "4.75 " << x4 << std::endl;
  std::cout << "DisY5:   " << "5.0  " << x5 << std::endl;
  std::cout << "DisY6:   " << "5.25 " << x6 << std::endl;
  std::cout << "DisY7:   " << "5.5  " << x7 << std::endl;
  std::cout << "DisY8:   " << "5.75 " << x8 << std::endl;
  std::cout << "DisY9:   " << "6.0  " << x9 << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << std::endl;

}



