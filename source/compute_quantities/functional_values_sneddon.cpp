#include "./../../include/dynamic_fracture.h"


template <int dim>
void
Dynamic_Fracture_Problem<dim>::compute_functional_values_Sneddon ()
{
  if (dim == 2)
    {

      double lines_sneddon_2d[] =
      {
	0.0, 1.0, 1.5, 
	1.78125, 1.8125, 1.875, 1.9375, 
	2.0,
	2.0625, 2.125, 2.1875, 2.21875, 
	2.5, 3.0, 4.0
      };
      // 2.15625

      const unsigned int n_lines_sneddon_2d = sizeof(lines_sneddon_2d) / sizeof(*lines_sneddon_2d);

      static unsigned int no = 0;
      ++no;
      std::ostringstream filename;
      filename << filename_basis_cod << Utilities::int_to_string(no, 2) << ".txt";
      std::cout << "writing " << filename.str() << std::endl;

      std::ofstream f(filename.str().c_str());

      if (test_case == "Sneddon")
        {
          for (unsigned int i = 0; i < n_lines_sneddon_2d; ++i)
            {
              double value = compute_cod(lines_sneddon_2d[i]);
              f << lines_sneddon_2d[i] << " " << value << std::endl;
            }
        }

    }
 
}
