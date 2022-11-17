#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
 
#include <iostream>
#include <fstream>
 
#include <map>
 
using namespace dealii;
 
 
template <int dim>
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string &       filename)
{
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;
 
  {
    std::map<types::boundary_id, unsigned int> boundary_count;
    for (const auto &face : triangulation.active_face_iterators())
      if (face->at_boundary())
        boundary_count[face->boundary_id()]++;
 
    std::cout << " boundary indicators: ";
    for (const std::pair<const types::boundary_id, unsigned int> &pair :
         boundary_count)
      {
        std::cout << pair.first << '(' << pair.second << " times) ";
      }
    std::cout << std::endl;
  }
 
  std::ofstream out(filename);
  GridOut       grid_out;
  grid_out.write_vtu(triangulation, out);
  std::cout << " written to " << filename << std::endl << std::endl;
}


void grid_1()
{

  Triangulation<2> triangulation;
 
  GridIn<2> gridin;
  gridin.attach_triangulation(triangulation);
  // std::string file_name = "mesh_files/mesh-new-1.msh";
  // std::string file_name = "mesh_files/mesh-1.msh";
  // std::string file_name = "mesh_files/gmsh/mesh.inp";
  std::string file_name = "mesh_files/gmsh/mesh.msh";
  std::ifstream f(file_name.c_str());
  if (file_name.find(".inp")!= std::string::npos)
    { 
      gridin.read_ucd(f);
      // gridin.read_abaqus(f);
      }
  else if (file_name.find(".msh")!= std::string::npos)
    gridin.read_msh(f);
  else
    AssertThrow(false, ExcMessage("Unknown file format"));
  // gridin.read_msh(f);
 
  
 
  print_mesh_info(triangulation, "grid-1.vtu");
}

int main()
{
  try
    {
      grid_1();

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
 
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}