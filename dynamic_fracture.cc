// Phase-field fracturetimer_newton
// Properties:
// - computes fracture propagation in elasticity
// - fully-coupled displacement phase-field model
// - predictor-corrector mesh adaptivity for small eps
// - tuncation-based goal-oriented adaptive time step control

// TODO: goal functionals
// augmented Lagrangian stopping only phase-field

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>  

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>

// #include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>


// C++
#include <fstream>
#include <sstream>

#include "include/utils.h"

// At the end of this top-matter, we import
// all deal.II names into the global
// namespace:				
using namespace dealii;
using namespace std;




// For Example 3 (multiple cracks in a heterogenous medium)
// reads .pgm file and returns it as floating point values
// taken from step-42
class BitmapFile
{
public:
  BitmapFile(const std::string &name);

  double
  get_value(const double x, const double y) const;

private:
  std::vector<double> image_data;
  double hx, hy;
  int nx, ny;

  double
  get_pixel_value(const int i, const int j) const;
};

// The constructor of this class reads in the data that describes
// the obstacle from the given file name.
BitmapFile::BitmapFile(const std::string &name)
  :
  image_data(0),
  hx(0),
  hy(0),
  nx(0),
  ny(0)
{
  std::ifstream f(name.c_str());
  AssertThrow (f, ExcMessage (std::string("Can't read from file <") +
                              name + ">!"));


  std::string temp;
  getline(f, temp);
  f >> temp;
  if (temp[0]=='#')
    getline(f, temp);

  f >> nx >> ny;

  AssertThrow(nx > 0 && ny > 0, ExcMessage("Invalid file format."));

  for (int k = 0; k < nx * ny; k++)
    {
      unsigned int val;
      f >> val;
      image_data.push_back(val / 255.0);
    }

  hx = 1.0 / (nx - 1);
  hy = 1.0 / (ny - 1);
}

// The following two functions return the value of a given pixel with
// coordinates $i,j$, which we identify with the values of a function
// defined at positions <code>i*hx, j*hy</code>, and at arbitrary
// coordinates $x,y$ where we do a bilinear interpolation between
// point values returned by the first of the two functions. In the
// second function, for each $x,y$, we first compute the (integer)
// location of the nearest pixel coordinate to the bottom left of
// $x,y$, and then compute the coordinates $\xi,\eta$ within this
// pixel. We truncate both kinds of variables from both below
// and above to avoid problems when evaluating the function outside
// of its defined range as may happen due to roundoff errors.
double
BitmapFile::get_pixel_value(const int i,
                            const int j) const
{
  assert(i >= 0 && i < nx);
  assert(j >= 0 && j < ny);
  return image_data[nx * (ny - 1 - j) + i];
}

double
BitmapFile::get_value(const double x,
                      const double y) const
{
  const int ix = std::min(std::max((int) (x / hx), 0), nx - 2);
  const int iy = std::min(std::max((int) (y / hy), 0), ny - 2);

  const double xi  = std::min(std::max((x-ix*hx)/hx, 1.), 0.);
  const double eta = std::min(std::max((y-iy*hy)/hy, 1.), 0.);

  return ((1-xi)*(1-eta)*get_pixel_value(ix,iy)
          +
          xi*(1-eta)*get_pixel_value(ix+1,iy)
          +
          (1-xi)*eta*get_pixel_value(ix,iy+1)
          +
          xi*eta*get_pixel_value(ix+1,iy+1));
}

template <int dim>
class BitmapFunction : public Function<dim>
{
public:
  BitmapFunction(const std::string &filename,
                 double x1_, double x2_, double y1_, double y2_, double minvalue_, double maxvalue_)
    : Function<dim>(1),
      f(filename), x1(x1_), x2(x2_), y1(y1_), y2(y2_), minvalue(minvalue_), maxvalue(maxvalue_)
  {}

  virtual
  double value (const Point<dim> &p,
                const unsigned int /*component*/ = 0) const
  {
    //Assert(dim==2, ExcNotImplemented());
    double x = (p(0)-x1)/(x2-x1);
    double y = (p(1)-y1)/(y2-y1);

    if (dim == 2)
      return minvalue + f.get_value(x,y)*(maxvalue-minvalue);
    else if (dim == 3)
      {
        double z = (p(2)-y1)/(y2-y1);
        return minvalue + (f.get_value(x,y)+f.get_value(x,z)+f.get_value(z,y))*(maxvalue-minvalue)/3.0;
      }

  }
private:
  BitmapFile f;
  double x1,x2,y1,y2;
  double minvalue, maxvalue;
};





 







///////////////////////////////////////////////////////////////////////////
// Class for initial values
  template <int dim>
  class InitialValuesPressurized : public Function<dim>
  {
    public:
      InitialValuesPressurized (const double alpha_eps) : Function<dim>(dim+dim+1) 
    {
      _alpha_eps = alpha_eps;
    }

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;

  private:
    double _alpha_eps;

  };


  template <int dim>
  double
  InitialValuesPressurized<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
    // only phase field
    double top = 2.0 +    _alpha_eps;
    double bottom = 2.0 - _alpha_eps;

    if (component == 2)   
      {
	//return 1.0; 
	if (((p(0) >= 1.5) && (p(0) <= 2.5)) &&
	    ((p(1) >= bottom) && (p(1) <= top))
	    )
	  {
	    return 0.0; 

	  }
	else 
	  return 1.0;
      }
    
    return 0.0;
  }


  template <int dim>
  void
  InitialValuesPressurized<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesPressurized<dim>::value (p, comp);
  }


///////////////////////////////////////////////////////////////////////////
// Class for initial values
  template <int dim>
  class InitialValuesSneddon3D : public Function<dim>
  {
    public:
      InitialValuesSneddon3D (const double alpha_eps) : Function<dim>(dim+1+dim) 
    {
      _alpha_eps = alpha_eps;
    }

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;

  private:
    double _alpha_eps;

  };


  template <int dim>
  double
  InitialValuesSneddon3D<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
 
    double top = 5.0 + _alpha_eps/2.0;
    double bottom = 5.0 - _alpha_eps/2.0;
    double radius = 1.0;
    // only phase field
    if (component == dim)   
      {
	//return 1.0; 
	if ((((p(0) - 5.0) * (p(0) - 5.0) + (p(2) - 5.0) * (p(2) - 5.0)) <=  radius * radius) &&
	    ((p(1) >= bottom)  && (p(1) <= top))
	    )
	  return 0.0;
	else 
	  return 1.0;
      }
    
    return 0.0;

  }


  template <int dim>
  void
  InitialValuesSneddon3D<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesSneddon3D<dim>::value (p, comp);
  }



///////////////////////////////////////////////////////////////////////////
// Class for initial values
  template <int dim>
  class InitialValuesHet3D : public Function<dim>
  {
    public:
      InitialValuesHet3D (const double alpha_eps) : Function<dim>(dim+1+dim) 
    {
      _alpha_eps = alpha_eps;
    }

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;

  private:
    double _alpha_eps;

  };


  template <int dim>
  double
  InitialValuesHet3D<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
 
    double width = _alpha_eps/2.0;
    double top = 5.0 + _alpha_eps/2.0;
    double bottom = 5.0 - _alpha_eps/2.0;
    double radius = 1.0;
    // only phase field
    if (component == dim)   
      {
	if (((p(0) >= 2.6 - width/2.0) && (p(0) <= 2.6 + width/2.0))
	    && ((p(1) >= 3.8- width/2.0) && (p(1) <= 5.5 + width/2.0))
	    && (p(2) >=4 - width/2.0) && (p(2) <=4 + width/2.0)
	    )
	  return 0.0;
	else if (((p(0) >= 5.5) && (p(0) <= 7.0))
		 && ((p(1) >= 4.0 - width/2.0) && (p(1) <= 4.0 + width/2.0))
		 && (p(2) >=6 - width/2.0) && (p(2) <=6 + width/2.0)
		 )
	  return 0.0;
	else
	  return 1.0;


      }
    
    return 0.0;

  }


  template <int dim>
  void
  InitialValuesHet3D<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesHet3D<dim>::value (p, comp);
  }





///////////////////////////////////////////////////////////////////////////
// Class for initial values
  template <int dim>
  class InitialValuesMiehe : public Function<dim>
  {
    public:
      InitialValuesMiehe (const double alpha_eps, 
			  const bool bool_initial_crack_via_phase_field) : Function<dim>(dim+1+dim) 
    {
      _alpha_eps = alpha_eps;
_bool_initial_crack_via_phase_field = bool_initial_crack_via_phase_field;
    }

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;

  private:
    double _alpha_eps;
    bool _bool_initial_crack_via_phase_field;

  };


  template <int dim>
  double
  InitialValuesMiehe<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
    // only phase field
	if (_bool_initial_crack_via_phase_field)
	{
    double top = 0.5 + _alpha_eps/2.0;
    double bottom = 0.5 - _alpha_eps/2.0;

    if (component == 2)   
      {
	//return 1.0; 
	if (((p(0) >= 0.5) && (p(0) <= 1.0)) &&
	    ((p(1) >= bottom) && (p(1) <= top))
	    )
	  {
	    return 0.0; 

	  }
	else 
	  return 1.0;
      }
}
	else
{
 if (component == 2)   
      {
	return 1.0;
      }

}

    
    return 0.0;


  }


  template <int dim>
  void
  InitialValuesMiehe<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesMiehe<dim>::value (p, comp);
  }




///////////////////////////////////////////////////////////////////////////
// Class for initial values
  template <int dim>
  class InitialValuesHohlwalzung : public Function<dim>
  {
    public:
      InitialValuesHohlwalzung (const double alpha_eps) : Function<dim>(dim+1+dim) 
    {
      _alpha_eps = alpha_eps;
    }

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;

  private:
    double _alpha_eps;

  };


  template <int dim>
  double
  InitialValuesHohlwalzung<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
  {
  
    if (component == dim)
      {
	if (((p(0) >= 0.0 - _alpha_eps) && (p(0) <= 0.0 + _alpha_eps)) &&
	    ((p(1) >= (-13.0)) && (p(1) <= (-7.0))) // -13.0 => 6mm
	    )
	  return 0.0;
	else 
	  return 1.0;
      }
    

    return 0.0;
  }


  template <int dim>
  void
  InitialValuesHohlwalzung<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesHohlwalzung<dim>::value (p, comp);
  }





////////////////////////////////////////////////////////////


template <int dim>
class NonhomDirichletBoundaryValues : public Function<dim> 
{
  public:
  NonhomDirichletBoundaryValues (const double time,
				 std::string test_case,
				 const double alpha_eps)    
    : Function<dim>(dim+1+dim) 
    {
      _time = time;  
      _test_case = test_case;
      _alpha_eps = alpha_eps;
    }
    
  virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const;

  virtual void vector_value (const Point<dim> &p, 
			     Vector<double>   &value) const;

private:
  double _time, _alpha_eps;
  std::string _test_case;

};

// The boundary values are given to component 
// with number 1.
template <int dim>
double
NonhomDirichletBoundaryValues<dim>::value (const Point<dim>  &p,
			     const unsigned int component) const
{
  Assert (component < this->n_components,
	  ExcIndexRange (component, 0, this->n_components));

  // TODO: scaled here and timestep * theta is commented
  double dis_step_per_timestep = 1.0; //1.0e-4;
  
  if (_test_case == "miehe_tension")
    {
      // Miehe tension
      if (component == 1) // u_y
	{
	  return ( ((p(1) == 1.0) && (p(0) <= 1.0) && (p(0) >= 0.0)) 
		   ? 			
		   (1.0) * _time * dis_step_per_timestep : 0 ); 
	  
	}
    }
  else if (_test_case == "miehe_shear")
    {
      // Miehe shear
      if (component == 0)
	{
	  return ( ((p(1) == 1.0) )
		   ?
		   (-1.0) * _time *dis_step_per_timestep : 0 );
	}

      if (component == dim+1)
	{
	  return ( ((p(1) == 1.0) )
		   ?
		   (-1.0) * dis_step_per_timestep : 0 );
	}


  // changed
  /*
      if (component == 2)
	{
	  double top = 0.5 + _alpha_eps/2.0;
	  double bottom = 0.5 - _alpha_eps/2.0;


	  if (p(0) == 1.0)
	    {
	      // Option 1 (sharp front)
//	      if ((p(1) >= bottom) && (p(1) <= top))
//		{
//		  return 0.0; 
//		}
//	      else 
//		return 1.0;

	      // Option 2 (smooth front of size eps)
	      return 1.0 + (-1.0) * std::exp(-std::abs(p(1) - 0.5)/(_alpha_eps));

	    }
	  
	} */

    } // bracket for mieher shear
    
 else if (_test_case == "l_shaped")
    {
      // L_shaped
      if (component == 1)
	{
	   
	  //  Geometry Jan 7, 2016
	  if ((p(1) == 250.0) && (p(0) <= 500.0) && (p(0) >= 470.0 ))
	    {
	      if (_time < 0.3)
		return _time * dis_step_per_timestep;
	      else if (_time >= 0.3 && _time < 0.8)
		return (0.6 - _time) * dis_step_per_timestep;
	      else if (_time >= 0.8)
		return (-1.0 + _time) * dis_step_per_timestep;


	    }
	  

	  /*
	  // Geometry Jan 8, 2016
	  if ((p(1) == 230.0) && (p(0) <= 470.0) && (p(0) >= 430.0 ))
	    {
	      if (_time < 0.3)
		return _time * dis_step_per_timestep;
	      else if (_time >= 0.3 && _time < 0.8)
		return (0.6 - _time) * dis_step_per_timestep;
	      else if (_time >= 0.8)
		return (-1.0 + _time) * dis_step_per_timestep;


	    }
	  */



	  // Old implementation
//	  return ( ((p(1) == 250.0) && (p(0) <= 500.0) && (p(0) >= 470.0 ))
//		   ?
//		   (1.0) * _time *dis_step_per_timestep : 0 );
	}


    } // end l-shaped
 else if (_test_case == "screw_domi")
   {
     if (component == 1)
       {
	 return ( ((p(1) == 0.0) )
		  ?
		  (1.0) * _time *dis_step_per_timestep : 0 );
       }
     //std::cout << "Bin drin." << std::endl;

   }



 
  return 0;
}



template <int dim>
void
NonhomDirichletBoundaryValues<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const 
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values (c) = NonhomDirichletBoundaryValues<dim>::value (p, c);
}


//------------------------------------------------------------//
 template <int dim>
  struct PointHistory
  {
    Tensor<2,dim> old_stress;
  };







/////////////////////////////////////////////////////////////////////
template <int dim>
class Dynamic_Fracture_Problem 
{
public:
  
  Dynamic_Fracture_Problem (const unsigned int degree);
  ~Dynamic_Fracture_Problem (); 
  void run ();
  
private:
  
  void set_global_parameters ();

  void set_runtime_parameters_Miehe ();
  void set_runtime_parameters_L_shaped ();
  void set_runtime_parameters_Sneddon ();
  void set_runtime_parameters_pressurized ();
  void set_runtime_parameters_screw_domi ();
  void set_runtime_parameters_Sneddon3D ();
  void set_runtime_parameters_Het3D ();

  void setup_system ();
  void make_material_vectors();
  void assemble_system_matrix ();   
  void assemble_system_rhs ();
  
  // functions for setting initial and boundary condition
  void set_initial_bc (const double time);
  void set_newton_bc ();
  
  void solve ();
  double newton_iteration(const double time);

  void solve_spatial_problem();

 // Nonlinear solver: error-based descend
  void newton_iteration_error_based(const double time);	

		  
  void output_results (const unsigned int refinement_cycle,
		       const BlockVector<double> solution) const;
  
  // Compute functional values
  double compute_point_value (Point<dim> p,
			      const unsigned int component) const;
  
  void compute_functional_values ();
 

  void compute_functional_values_Sneddon (); 
  double compute_cod (const double eval_line); 
  void compute_cod_Sneddon3D();

  void compute_stress_per_cell ();

  void compute_energy ();

  double goal_functional_stress_x ();

  // Mesh refinement
  bool refine_mesh();
  void project_back_phase_field (); 


  void setup_quadrature_point_history ();
  void update_quadrature_point_history ();

  std::vector<PointHistory<dim> > quadrature_point_history;

  const unsigned int   degree;
  
  Triangulation<dim>   triangulation;
  FESystem<dim>        fe;
  DoFHandler<dim>      dof_handler;

  ConstraintMatrix     constraints;
  
  BlockSparsityPattern      sparsity_pattern; 
  BlockSparseMatrix<double> system_matrix; 
  
  BlockVector<double> solution, newton_update, old_timestep_solution, 
    old_old_timestep_solution, tmp_solution_for_time_adaptivity;
  BlockVector<double> system_rhs;
  BlockVector<double> solution_lambda_penal_func, old_timestep_solution_lambda_penal_func;


  Vector<double> lame_coefficient_mu_vector, lame_coefficient_lambda_vector;
  Vector<float> solution_stress_per_cell;

  TimerOutput         timer;
  // ConditionalOStream pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_com) == 0));

  // MPI_Comm mpi_com;

  Function<dim> *func_emodulus;
  
  // Global variables for timestepping scheme   
  unsigned int timestep_number;
  unsigned int max_no_timesteps;  
  double timestep, theta, time; 
  std::string time_stepping_scheme;
  std::string test_case, sub_test_case;
  double old_timestep, old_old_timestep;

  double force_structure_x, force_structure_y;	  
  

  double gravity_x, gravity_y, volume_source, traction_x, traction_y;
  
  
  // Structure parameters
  double density_structure; 
  double lame_coefficient_mu, lame_coefficient_lambda, poisson_ratio_nu;  


  double cell_diameter;  
  double max_no_of_augmented_L_penal_iterations;
  unsigned int penal_iterations;
  
  bool bool_use_stress_splitting;
 
  double constant_k, alpha_eps, G_c, delta_penal, gamma_penal;
  double upper_newton_rho;
  double tolerance_augmented_L_iteration,  tolerance_absolute_augmented_L_iteration;
  double current_pressure, old_timestep_current_pressure, alpha_biot;
  double E_modulus, E_prime;
  unsigned int total_number_newton_steps;


  double lower_bound_newton_update;
  double lower_bound_newton_residuum; 

  bool bool_use_pf_extra, bool_use_error_oriented_Newton, bool_use_modified_Newton, 
    bool_use_adaptive_newton_bound, bool_plot_additional_solutions,
    bool_set_explicitely_delta_fp, bool_set_explicitely_constant_k;
  unsigned int max_no_newton_steps, max_no_line_search_steps;
  std::string filename_basis, filename_basis_cod;

  double value_phase_field_for_refinement;

  unsigned int global_refinement_steps;
  unsigned int pred_corr_levels;

  double TOL_for_time_adaptivity, timestep_rejection_factor, timestep_growth_factor;
  bool use_time_step_control;
  double timestep_upper_bound, timestep_lower_bound,end_time_value;
  int number_of_nonadaptive_time_steps;

  double delta_fixed_point_newton;
  double a_fp, b_fp;

  bool bool_set_initial_strain_history, bool_initial_crack_via_phase_field, bool_use_strain_history;

  bool bool_use_dynamic_code_with_velocities;
  

};


// The constructor of this class is comparable 
// to other tutorials steps, e.g., step-22, and step-31.
template <int dim>
Dynamic_Fracture_Problem<dim>::Dynamic_Fracture_Problem (const unsigned int degree)
                :
                degree (degree),
		triangulation (Triangulation<dim>::maximum_smoothing),
                fe (FE_Q<dim>(degree), dim,  // velocities                
		    FE_Q<dim>(degree), 1,    // phase-field
		    FE_Q<dim>(degree), dim),   // displacements
                dof_handler (triangulation),
    timer (std::cout, TimerOutput::summary, TimerOutput::cpu_times)
		// timer (std::cout, TimerOutput::summary, TimerOutput::cpu_times)
    // timer (pcout, TimerOutput::summary, TimerOutput::cpu_and_wall_times)
    // timer(mpi_com, pcout, TimerOutput::every_call_and_summary,
    //     TimerOutput::cpu_and_wall_times)		
{}


// This is the standard destructor.
template <int dim>
Dynamic_Fracture_Problem<dim>::~Dynamic_Fracture_Problem () 
{}


// Now, there follow several functions to perform
// the spectral decomposition of the stress tensor
// into tension and compression parts
// assumes the matrix is symmetric!
// The explicit calculation does only work
// in 2d. For 3d, we should use other libraries or approximative
// tools to compute eigenvectors and -functions.
// Borden et al. (2012, 2013) suggested some papers to look into.
template <int dim>
void eigen_vectors_and_values(
  double &E_eigenvalue_1, double &E_eigenvalue_2,
  Tensor<2,dim> &ev_matrix,
  const Tensor<2,dim> &matrix)
{
  double sq = std::sqrt((matrix[0][0] - matrix[1][1]) * (matrix[0][0] - matrix[1][1]) + 4.0*matrix[0][1]*matrix[1][0]);
  E_eigenvalue_1 = 0.5 * ((matrix[0][0] + matrix[1][1]) + sq);
  E_eigenvalue_2 = 0.5 * ((matrix[0][0] + matrix[1][1]) - sq);

  // Compute eigenvectors
  Tensor<1,dim> E_eigenvector_1;
  Tensor<1,dim> E_eigenvector_2;
  if (std::abs(matrix[0][1]) < 1e-10*std::abs(matrix[0][0]))
    {
      // E is close to diagonal
      E_eigenvector_1[0]=0;
      E_eigenvector_1[1]=1;
      E_eigenvector_2[0]=1;
      E_eigenvector_2[1]=0;
    }
  else
    {
      E_eigenvector_1[0] = 1.0/(std::sqrt(1 + (E_eigenvalue_1 - matrix[0][0])/matrix[0][1] * (E_eigenvalue_1 - matrix[0][0])/matrix[0][1]));
      E_eigenvector_1[1] = (E_eigenvalue_1 - matrix[0][0])/(matrix[0][1] * (std::sqrt(1 + (E_eigenvalue_1 - matrix[0][0])/matrix[0][1] * (E_eigenvalue_1 - matrix[0][0])/matrix[0][1])));
      E_eigenvector_2[0] = 1.0/(std::sqrt(1 + (E_eigenvalue_2 - matrix[0][0])/matrix[0][1] * (E_eigenvalue_2 - matrix[0][0])/matrix[0][1]));
      E_eigenvector_2[1] = (E_eigenvalue_2 - matrix[0][0])/(matrix[0][1] * (std::sqrt(1 + (E_eigenvalue_2 - matrix[0][0])/matrix[0][1] * (E_eigenvalue_2 - matrix[0][0])/matrix[0][1])));
    }

  ev_matrix[0][0] = E_eigenvector_1[0];
  ev_matrix[0][1] = E_eigenvector_2[0];
  ev_matrix[1][0] = E_eigenvector_1[1];
  ev_matrix[1][1] = E_eigenvector_2[1];

  // Sanity check if orthogonal
  double scalar_prod = 1.0e+10;
  scalar_prod = E_eigenvector_1[0] * E_eigenvector_2[0] + E_eigenvector_1[1] * E_eigenvector_2[1];

  if (scalar_prod > 1.0e-6)
    {
      // TODO: check this condition and also with the old code
      std::cout << "Seems not to be orthogonal" << std::endl;
      //abort();
    }

}

template <int dim>
void decompose_stress(
  Tensor<2,dim> &stress_term_plus,
  Tensor<2,dim> &stress_term_minus,
  const Tensor<2, dim> &E,
  const double tr_E,
  const Tensor<2, dim> &E_LinU,
  const double tr_E_LinU,
  const double lame_coefficient_lambda,
  const double lame_coefficient_mu,
  const bool derivative)
{
  // This is the identity matrix in two dimensions:
  const Tensor<2,dim> Identity = ALE_Transformations
    ::get_Identity<dim> ();

  Tensor<2, dim> zero_matrix;
  zero_matrix.clear();


  // Compute first the eigenvalues for u (as in the previous function)
  // and then for \delta u

  // Compute eigenvalues/vectors
  double E_eigenvalue_1, E_eigenvalue_2;
  Tensor<2,dim> P_matrix;
  eigen_vectors_and_values(E_eigenvalue_1, E_eigenvalue_2,P_matrix,E);

  double E_eigenvalue_1_plus = std::max(0.0, E_eigenvalue_1);
  double E_eigenvalue_2_plus = std::max(0.0, E_eigenvalue_2);

  Tensor<2,dim> Lambda_plus;
  Lambda_plus[0][0] = E_eigenvalue_1_plus;
  Lambda_plus[0][1] = 0.0;
  Lambda_plus[1][0] = 0.0;
  Lambda_plus[1][1] = E_eigenvalue_2_plus;

  if (!derivative)
    {
      Tensor<2,dim> E_plus = P_matrix * Lambda_plus * transpose(P_matrix);

      double tr_E_positive = std::max(0.0, tr_E);

      stress_term_plus = lame_coefficient_lambda * tr_E_positive * Identity
                         + 2 * lame_coefficient_mu * E_plus;

      stress_term_minus = lame_coefficient_lambda * (tr_E - tr_E_positive) * Identity
                          + 2 * lame_coefficient_mu * (E - E_plus);
    }
  else
    {
      // Derviatives (\delta u)

      // Compute eigenvalues/vectors
      double E_eigenvalue_1_LinU, E_eigenvalue_2_LinU;
      Tensor<1,dim> E_eigenvector_1_LinU;
      Tensor<1,dim> E_eigenvector_2_LinU;
      Tensor<2,dim> P_matrix_LinU;

      // Compute linearized Eigenvalues
      double diskriminante = std::sqrt(E[0][1] * E[1][0] + (E[0][0] - E[1][1]) * (E[0][0] - E[1][1])/4.0);

      E_eigenvalue_1_LinU = 0.5 * tr_E_LinU + 1.0/(2.0 * diskriminante) *
                            (E_LinU[0][1] * E[1][0] + E[0][1] * E_LinU[1][0] + (E[0][0] - E[1][1])*(E_LinU[0][0] - E_LinU[1][1])/2.0);

      E_eigenvalue_2_LinU = 0.5 * tr_E_LinU - 1.0/(2.0 * diskriminante) *
                            (E_LinU[0][1] * E[1][0] + E[0][1] * E_LinU[1][0] + (E[0][0] - E[1][1])*(E_LinU[0][0] - E_LinU[1][1])/2.0);


      // Compute normalized Eigenvectors and P
      double normalization_1 = 1.0/(std::sqrt(1 + (E_eigenvalue_1 - E[0][0])/E[0][1] * (E_eigenvalue_1 - E[0][0])/E[0][1]));
      double normalization_2 = 1.0/(std::sqrt(1 + (E_eigenvalue_2 - E[0][0])/E[0][1] * (E_eigenvalue_2 - E[0][0])/E[0][1]));

      double normalization_1_LinU = 0.0;
      double normalization_2_LinU = 0.0;

      normalization_1_LinU = -1.0 * (1.0/(1.0 + (E_eigenvalue_1 - E[0][0])/E[0][1] * (E_eigenvalue_1 - E[0][0])/E[0][1])
                                     * 1.0/(2.0 * std::sqrt(1.0 + (E_eigenvalue_1 - E[0][0])/E[0][1] * (E_eigenvalue_1 - E[0][0])/E[0][1]))
                                     * (2.0 * (E_eigenvalue_1 - E[0][0])/E[0][1])
                                     * ((E_eigenvalue_1_LinU - E_LinU[0][0]) * E[0][1] - (E_eigenvalue_1 - E[0][0]) * E_LinU[0][1])/(E[0][1] * E[0][1]));

      normalization_2_LinU = -1.0 * (1.0/(1.0 + (E_eigenvalue_2 - E[0][0])/E[0][1] * (E_eigenvalue_2 - E[0][0])/E[0][1])
                                     * 1.0/(2.0 * std::sqrt(1.0 + (E_eigenvalue_2 - E[0][0])/E[0][1] * (E_eigenvalue_2 - E[0][0])/E[0][1]))
                                     * (2.0 * (E_eigenvalue_2 - E[0][0])/E[0][1])
                                     * ((E_eigenvalue_2_LinU - E_LinU[0][0]) * E[0][1] - (E_eigenvalue_2 - E[0][0]) * E_LinU[0][1])/(E[0][1] * E[0][1]));


      E_eigenvector_1_LinU[0] = normalization_1 * 1.0;
      E_eigenvector_1_LinU[1] = normalization_1 * (E_eigenvalue_1 - E[0][0])/E[0][1];

      E_eigenvector_2_LinU[0] = normalization_2 * 1.0;
      E_eigenvector_2_LinU[1] = normalization_2 * (E_eigenvalue_2 - E[0][0])/E[0][1];


      // Apply product rule to normalization and vector entries
      double EV_1_part_1_comp_1 = 0.0;  // LinU in vector entries, normalization U
      double EV_1_part_1_comp_2 = 0.0;  // LinU in vector entries, normalization U
      double EV_1_part_2_comp_1 = 0.0;  // vector entries U, normalization LinU
      double EV_1_part_2_comp_2 = 0.0;  // vector entries U, normalization LinU

      double EV_2_part_1_comp_1 = 0.0;  // LinU in vector entries, normalization U
      double EV_2_part_1_comp_2 = 0.0;  // LinU in vector entries, normalization U
      double EV_2_part_2_comp_1 = 0.0;  // vector entries U, normalization LinU
      double EV_2_part_2_comp_2 = 0.0;  // vector entries U, normalization LinU

      // Effizienter spaeter, aber erst einmal uebersichtlich und verstehen!
      EV_1_part_1_comp_1 = normalization_1 * 0.0;
      EV_1_part_1_comp_2 = normalization_1 *
                           ((E_eigenvalue_1_LinU - E_LinU[0][0]) * E[0][1] - (E_eigenvalue_1 - E[0][0]) * E_LinU[0][1])/(E[0][1] * E[0][1]);

      EV_1_part_2_comp_1 = normalization_1_LinU * 1.0;
      EV_1_part_2_comp_2 = normalization_1_LinU * (E_eigenvalue_1 - E[0][0])/E[0][1];


      EV_2_part_1_comp_1 = normalization_2 * 0.0;
      EV_2_part_1_comp_2 = normalization_2 *
                           ((E_eigenvalue_2_LinU - E_LinU[0][0]) * E[0][1] - (E_eigenvalue_2 - E[0][0]) * E_LinU[0][1])/(E[0][1] * E[0][1]);

      EV_2_part_2_comp_1 = normalization_2_LinU * 1.0;
      EV_2_part_2_comp_2 = normalization_2_LinU * (E_eigenvalue_2 - E[0][0])/E[0][1];



      // Build eigenvectors
      E_eigenvector_1_LinU[0] = EV_1_part_1_comp_1 + EV_1_part_2_comp_1;
      E_eigenvector_1_LinU[1] = EV_1_part_1_comp_2 + EV_1_part_2_comp_2;

      E_eigenvector_2_LinU[0] = EV_2_part_1_comp_1 + EV_2_part_2_comp_1;
      E_eigenvector_2_LinU[1] = EV_2_part_1_comp_2 + EV_2_part_2_comp_2;



      // P-Matrix
      P_matrix_LinU[0][0] = E_eigenvector_1_LinU[0];
      P_matrix_LinU[0][1] = E_eigenvector_2_LinU[0];
      P_matrix_LinU[1][0] = E_eigenvector_1_LinU[1];
      P_matrix_LinU[1][1] = E_eigenvector_2_LinU[1];


      double E_eigenvalue_1_plus_LinU = 0.0;
      double E_eigenvalue_2_plus_LinU = 0.0;


      // Very important: Set E_eigenvalue_1_plus_LinU to zero when
      // the corresponding rhs-value is set to zero and NOT when
      // the value itself is negative!!!
      if (E_eigenvalue_1 < 0.0)
        {
          E_eigenvalue_1_plus_LinU = 0.0;
        }
      else
        E_eigenvalue_1_plus_LinU = E_eigenvalue_1_LinU;


      if (E_eigenvalue_2 < 0.0)
        {
          E_eigenvalue_2_plus_LinU = 0.0;
        }
      else
        E_eigenvalue_2_plus_LinU = E_eigenvalue_2_LinU;



      Tensor<2,dim> Lambda_plus_LinU;
      Lambda_plus_LinU[0][0] = E_eigenvalue_1_plus_LinU;
      Lambda_plus_LinU[0][1] = 0.0;
      Lambda_plus_LinU[1][0] = 0.0;
      Lambda_plus_LinU[1][1] = E_eigenvalue_2_plus_LinU;

      Tensor<2,dim> E_plus_LinU = P_matrix_LinU * Lambda_plus * transpose(P_matrix) +  P_matrix * Lambda_plus_LinU * transpose(P_matrix) + P_matrix * Lambda_plus * transpose(P_matrix_LinU);


      double tr_E_positive_LinU = 0.0;
      if (tr_E < 0.0)
        {
          tr_E_positive_LinU = 0.0;

        }
      else
        tr_E_positive_LinU = tr_E_LinU;



      stress_term_plus = lame_coefficient_lambda * tr_E_positive_LinU * Identity
                         + 2 * lame_coefficient_mu * E_plus_LinU;

      stress_term_minus = lame_coefficient_lambda * (tr_E_LinU - tr_E_positive_LinU) * Identity
                          + 2 * lame_coefficient_mu * (E_LinU - E_plus_LinU);


      // Sanity check
      //Tensor<2,dim> stress_term = lame_coefficient_lambda * tr_E_LinU * Identity
      //  + 2 * lame_coefficient_mu * E_LinU;

      //std::cout << stress_term.norm() << "   " << stress_term_plus.norm() << "   " << stress_term_minus.norm() << std::endl;
    }


}

// In this method, we set up runtime parameters 
template <int dim>
void Dynamic_Fracture_Problem<dim>::set_global_parameters ()
{
  // Initialization of some variables
  // that may be overwritten in the set_runtime routines
  bool_set_explicitely_delta_fp = false;
  bool_set_explicitely_constant_k = false;
  
  bool_use_strain_history = false;
  bool_initial_crack_via_phase_field = false;

  upper_newton_rho = 0.9;

  constant_k = 1.1e-2; // it's kappa
  alpha_eps = 2.2e-2; // it's eps

  use_time_step_control   = false;

  bool_use_dynamic_code_with_velocities = true;

}


// In this method, we set up runtime parameters 
template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_Miehe ()
{
  // it does not have anything with boundary condition
  // changed bool_initial_crack_via_phase_field = false;
  bool_initial_crack_via_phase_field = false;
  

  // Parameters
  current_pressure = 0.0; 
  alpha_biot = 0.0;

  G_c = 2.7;
  delta_penal = 0.0; // simple penalization
  //changed gamma_penal =  1.0 
  gamma_penal = 1.0e+3; //  augmented Lagrangian penalization

  gravity_x = 0.0;
  gravity_y = 0.0;
  volume_source = 0.0;


  // Elasticity parameters
  force_structure_x = 0.0;
  force_structure_y = 0.0;

  // Traction 
  traction_x = 0.0;
  traction_y = 0.0; 


  density_structure = 1.0; 

  // Structure parameters
  lame_coefficient_mu = 80.77e+3; 
  poisson_ratio_nu = 0.3; 
  
  lame_coefficient_lambda = 121.15e+3; //(2 * poisson_ratio_nu * lame_coefficient_mu)/(1.0 - 2 * poisson_ratio_nu);

  
  // Timestepping schemes
  //BE, CN, CN_shifted
  // BE: backward euler
  time_stepping_scheme = "BE";

  // Timestep size:
  timestep = 1.0e-5; //1.0e-4;

  // Maximum number of timesteps:
  max_no_timesteps = 1000; //130;
  // changed end_time_value = 1.0e-2;
  end_time_value = 1.3e-2; // Crack reaches lower left around 1.3e-2 sec

  number_of_nonadaptive_time_steps = 2;

  TOL_for_time_adaptivity = 1.0e-2; // 1.0e-2
  use_time_step_control   = false;


 
  timestep_rejection_factor = 0.5; // 0.1, 0.2 (Turek's suggestion)
  timestep_growth_factor = 1.0e+4; // 2.0
 
  timestep_upper_bound = 5.0e-5;
  timestep_lower_bound = 5.0e-6;



  
  // A variable to count the number of time steps
  timestep_number = 0;

  // Counts total time  
  time = 0;
 
  // Here, we choose a time-stepping scheme that
  // is based on finite differences:
  // BE         = backward Euler scheme 
  // CN         = Crank-Nicolson scheme
  // CN_shifted = time-shifted Crank-Nicolson scheme 
  // For further properties of these schemes,
  // we refer to standard literature.
  if (time_stepping_scheme == "BE")
    theta = 1.0;
  else if (time_stepping_scheme == "CN")
    theta = 0.5;
  else if (time_stepping_scheme == "CN_shifted")
    theta = 0.5 + timestep;
  else 
    std::cout << "No such timestepping scheme" << std::endl;

  // In the following, we read a *.inp grid from a file.
  std::string grid_name;
  if (!bool_initial_crack_via_phase_field)
    grid_name  = "unit_slit.inp";
  else if (bool_initial_crack_via_phase_field)
    grid_name  = "unit_square_1.inp"; 
 
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==2, ExcInternalError());
  grid_in.read_ucd (input_file); 
  
  global_refinement_steps = 4;
  pred_corr_levels = 0;
  triangulation.refine_global (global_refinement_steps); 

  //filename_basis  = "solution_Miehe_eps_2h_ref_6_delta_0_"; 
  filename_basis  = "solution_Miehe_test_"; 
  bool_use_error_oriented_Newton = false;
  bool_use_modified_Newton = true; // if true need to set error_oriented_Newton to false
  bool_set_explicitely_delta_fp = false; // if true, must set use_modified Newton to false
  bool_set_explicitely_constant_k = false;
  bool_use_pf_extra = false;

  a_fp = 1e-2; //0.001; //0.05;
  b_fp = 5.0; //1.5;
  max_no_of_augmented_L_penal_iterations = 10; // 20
  if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

  // changed tolerance_augmented_L_iteration = 1.0e-5;
  tolerance_augmented_L_iteration = 1.0e-2;
  tolerance_absolute_augmented_L_iteration = 1.0e-5;

  max_no_newton_steps  = 100;

  // Line search
  max_no_line_search_steps = 10; 
  if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  // Redefined below !!
  lower_bound_newton_update = 1.0e-10;
  lower_bound_newton_residuum = 1.0e-10; 

  // When `true' be careful because tolerances are fixed 
  //  and differ from the choices here.
  bool_use_adaptive_newton_bound = false;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();

  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = true;

}


// In this method, we set up runtime parameters 
template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_L_shaped ()
{

  // Parameters
  current_pressure = 0.0; 
  alpha_biot = 0.0;

  G_c = 8.9e-2;
  delta_penal = 0.0; // simple penalization
  gamma_penal = 1.0e-3; // augmented Lagrangian penalization

  gravity_x = 0.0;
  gravity_y = 0.0;
  volume_source = 0.0;


  // Elasticity parameters
  force_structure_x = 0.0;
  force_structure_y = 0.0;

  // Traction 
  traction_x = 0.0;
  traction_y = 0.0; 


  density_structure = 1.0; 

  // Structure parameters
  lame_coefficient_mu = 10.95e+3; 
  poisson_ratio_nu = 0.3; 
  
  lame_coefficient_lambda = 6.16e+3; //(2 * poisson_ratio_nu * lame_coefficient_mu)/(1.0 - 2 * poisson_ratio_nu);

  
  // Timestepping schemes
  //BE, CN, CN_shifted
  time_stepping_scheme = "BE";

  // Timestep size:
  timestep = 1.0e-3;

  // Maximum number of timesteps:
  max_no_timesteps = 2000;
  end_time_value = 1e+10;
 
  // A variable to count the number of time steps
  timestep_number = 0;

  // Counts total time  
  time = 0;
 
  // Here, we choose a time-stepping scheme that
  // is based on finite differences:
  // BE         = backward Euler scheme 
  // CN         = Crank-Nicolson scheme
  // CN_shifted = time-shifted Crank-Nicolson scheme 
  // For further properties of these schemes,
  // we refer to standard literature.
  if (time_stepping_scheme == "BE")
    theta = 1.0;
  else if (time_stepping_scheme == "CN")
    theta = 0.5;
  else if (time_stepping_scheme == "CN_shifted")
    theta = 0.5 + timestep;
  else 
    std::cout << "No such timestepping scheme" << std::endl;

  // In the following, we read a *.inp grid from a file.
  std::string grid_name;
  grid_name  = "l_shape_Jan_7_2016.inp"; 
  //grid_name  = "l_shape_Jan_8_2016.inp"; 
 
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==2, ExcInternalError());
  grid_in.read_ucd (input_file); 
  
  // TODO: insert all time step tolerances etc from Miehe oben
  // TODO: adjust goal functional
 
  global_refinement_steps = 2;
  pred_corr_levels = 0;   
  triangulation.refine_global (global_refinement_steps); 

  filename_basis  = "solution_l_shape_"; 
  bool_use_error_oriented_Newton = true;
  bool_use_pf_extra = true;

  max_no_of_augmented_L_penal_iterations = 20; // 50
 if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

  tolerance_augmented_L_iteration = 1.0e-4;
  tolerance_absolute_augmented_L_iteration = 1.0e-8;

  max_no_newton_steps  = 20; //50;
  max_no_line_search_steps = 10;
  if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  // reset in run method because of adaptive Newton technique
  lower_bound_newton_update = 1.0e-8;
  lower_bound_newton_residuum = 1.0e-8; 

  bool_use_adaptive_newton_bound = false;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();
 
  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = true;

}




// In this method, we set up runtime parameters 
template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_Sneddon ()
{

  // Parameters
  current_pressure = 1.0e-3; 
  alpha_biot = 0.0;

  G_c = 1.0;
  delta_penal = 0.0; // simple penalization
  gamma_penal = 1.0e+3; //1.0e+3; //1.0e+4; // augmented Lagrangian penalization

  gravity_x = 0.0;
  gravity_y = 0.0;
  volume_source = 0.0;


  // Elasticity parameters
  force_structure_x = 0.0;
  force_structure_y = 0.0;

  // Traction 
  traction_x = 0.0;
  traction_y = 0.0; 


  density_structure = 1.0; 

  // Structure parameters
  density_structure = 1.0; 
  poisson_ratio_nu = 0.2; 
  E_modulus = 1.0;
  E_prime = E_modulus/(1.0 - poisson_ratio_nu * poisson_ratio_nu);
  
  lame_coefficient_mu = E_modulus/(2.0*(1 + poisson_ratio_nu));
  
  lame_coefficient_lambda = (2 * poisson_ratio_nu * lame_coefficient_mu)/(1.0 - 2 * poisson_ratio_nu);

  
  // Timestepping schemes
  //BE, CN, CN_shifted
  time_stepping_scheme = "BE";

  // Timestep size:
  timestep = 1.0; //2.5e-1; //1.0;

  // Maximum number of timesteps:
  max_no_timesteps = 5;
  end_time_value = 1e+10;
  
  // A variable to count the number of time steps
  timestep_number = 0;

  // Counts total time  
  time = 0;
 
  // Here, we choose a time-stepping scheme that
  // is based on finite differences:
  // BE         = backward Euler scheme 
  // CN         = Crank-Nicolson scheme
  // CN_shifted = time-shifted Crank-Nicolson scheme 
  // For further properties of these schemes,
  // we refer to standard literature.
  if (time_stepping_scheme == "BE")
    theta = 1.0;
  else if (time_stepping_scheme == "CN")
    theta = 0.5;
  else if (time_stepping_scheme == "CN_shifted")
    theta = 0.5 + timestep;
  else 
    std::cout << "No such timestepping scheme" << std::endl;

  // In the following, we read a *.inp grid from a file.
  std::string grid_name;
  //grid_name  = "unit_slit.inp"; 
  // Example 2
  grid_name  = "unit_square_4.inp"; 
  
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==2, ExcInternalError());
  grid_in.read_ucd (input_file); 


 // TODO: insert all time step tolerances etc.
  
  global_refinement_steps = 6;
  pred_corr_levels = 2;
  triangulation.refine_global (global_refinement_steps); 

  //filename_basis  = "solution_pressurized_pf_extra_eps_2h_ref_7_adaptive_"; 
  filename_basis  = "solution_Sneddon_pf_extra_2h_ref_5_"; 
  filename_basis_cod = "cod_ref_6_";
  bool_use_error_oriented_Newton = false;
  bool_use_pf_extra = true;
  max_no_of_augmented_L_penal_iterations = 50; //15; //40;
 if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

  tolerance_augmented_L_iteration = 1.0e-3; // 1e-4;
  tolerance_absolute_augmented_L_iteration = 1.0e-8;

  max_no_newton_steps  = 20; //20; // 50
  max_no_line_search_steps = 10; 
  if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  lower_bound_newton_update = 1.0e-8;
  lower_bound_newton_residuum = 1.0e-8; 

  bool_use_adaptive_newton_bound = true;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();

  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = false;


}






// In this method, we set up runtime parameters 
template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_pressurized ()
{
  // Gotop // short cut for Thomas to go to this function

  // Info: double-check boundary conditions because
  // the trial and test functions change in the mixed dynamic
  // system in comparison to the quasi-static code
  // For the pressurized code bc should be fine since everything
  // is zero displacements. But for examples with non-homogeneous
  // displacements and velocities, this might be important
  // Change: step 2
  // bool_use_dynamic_code_with_velocities = true;
  bool_use_dynamic_code_with_velocities = false;

  // Parameters
  current_pressure = 1.0e-3; 
  alpha_biot = 0.0;

  G_c = 1.0;
  delta_penal = 0.0; // simple penalization
  gamma_penal = 1.0e+3; //1.0e+3; //1.0e+4; // augmented Lagrangian penalization

  gravity_x = 0.0;
  gravity_y = 0.0;
  volume_source = 0.0;


  // Elasticity parameters
  force_structure_x = 0.0;
  force_structure_y = 0.0;

  // Traction 
  traction_x = 0.0;
  traction_y = 0.0; 


  density_structure = 1.0; 

  // Structure parameters
  density_structure = 1.0; 
  poisson_ratio_nu = 0.2; 
  E_modulus = 1.0;
  E_prime = E_modulus/(1.0 - poisson_ratio_nu * poisson_ratio_nu);
  
  lame_coefficient_mu = E_modulus/(2.0*(1 + poisson_ratio_nu));
  
  lame_coefficient_lambda = (2 * poisson_ratio_nu * lame_coefficient_mu)/(1.0 - 2 * poisson_ratio_nu);

  
  // Timestepping schemes
  //BE, CN, CN_shifted
  time_stepping_scheme = "BE";

  // Timestep size:
  timestep = 1.0; //2.5e-1; //1.0;

  // Maximum number of timesteps:
  max_no_timesteps = 20;
  end_time_value = 1e+10;
  
  // A variable to count the number of time steps
  timestep_number = 0;

  // Counts total time  
  time = 0;
 
  // Here, we choose a time-stepping scheme that
  // is based on finite differences:
  // BE         = backward Euler scheme 
  // CN         = Crank-Nicolson scheme
  // CN_shifted = time-shifted Crank-Nicolson scheme 
  // For further properties of these schemes,
  // we refer to standard literature.
  if (time_stepping_scheme == "BE")
    theta = 1.0;
  else if (time_stepping_scheme == "CN")
    theta = 0.5;
  else if (time_stepping_scheme == "CN_shifted")
    theta = 0.5 + timestep;
  else 
    std::cout << "No such timestepping scheme" << std::endl;

  // In the following, we read a *.inp grid from a file.
  std::string grid_name;
  //grid_name  = "unit_slit.inp"; 
  // Example 2
  grid_name  = "unit_square_4.inp"; 
  
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==2, ExcInternalError());
  grid_in.read_ucd (input_file); 
  

 // TODO: insert all time step tolerances etc.
 
  global_refinement_steps = 6;
  pred_corr_levels = 0;   
  triangulation.refine_global (global_refinement_steps); 

  //filename_basis  = "solution_pressurized_pf_extra_eps_2h_ref_7_adaptive_"; 
  // filename_basis  = "solution_pressurized_dynamic_";  
  filename_basis  = "solution_pressurized_quasi_static_";
  bool_use_error_oriented_Newton = false;
  bool_use_modified_Newton = false; // if true and be used, then must set error oriented false
  bool_set_explicitely_delta_fp = false; // if true, must set use_modified Newton to false
  bool_set_explicitely_constant_k = false;
  bool_use_pf_extra = true;


  a_fp = 0.1; //0.05;
  b_fp = 2.0; //1.5;
  max_no_of_augmented_L_penal_iterations = 0; //15; //40;
 if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

  tolerance_augmented_L_iteration = 1.0e-2; // 1e-4;
  tolerance_absolute_augmented_L_iteration = 1.0e-8;

  max_no_newton_steps  = 100; //20; // 50
  max_no_line_search_steps = 10;
  if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  lower_bound_newton_update = 1.0e-8;
  lower_bound_newton_residuum = 1.0e-8; 

  bool_use_adaptive_newton_bound = false;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();

  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = false;
 
}


// In this method, we set up runtime parameters 
template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_screw_domi ()
{
  // screw
  //sub_test_case = "hohlwalzung";
 


  // Parameters
  current_pressure = 0.0; 
  alpha_biot = 0.0;

  G_c = 2.7;
  delta_penal = 0.0; // simple penalization
  gamma_penal = 1.0; // augmented Lagrangian penalization

  gravity_x = 0.0;
  gravity_y = 0.0;
  volume_source = 0.0;


  // Elasticity parameters
  force_structure_x = 0.0;
  force_structure_y = 0.0;

  // Traction 
  traction_x = 0.0;
  traction_y = 0.0; 


  density_structure = 1.0; 

  // Solid parameters
  lame_coefficient_mu = 80.77e+3;
  lame_coefficient_lambda = 121.15e+3;


  
  // Timestepping schemes
  //BE, CN, CN_shifted
  time_stepping_scheme = "BE";

  // Timestep size:
  timestep = 1.0e-2; //1.0e-4;

  // Maximum number of timesteps:
  max_no_timesteps = 20; //130;
  end_time_value = 1.0; // Crack reaches lower left around 1.3e-2 sec

  number_of_nonadaptive_time_steps = 2;

  TOL_for_time_adaptivity = 1.0e-2; // 1.0e-2
  use_time_step_control   = false;


 
  timestep_rejection_factor = 0.5; // 0.1, 0.2 (Turek's suggestion)
  timestep_growth_factor = 1.0e+4; // 2.0
 
  timestep_upper_bound = 5.0e-5;
  timestep_lower_bound = 5.0e-6;



  
  // A variable to count the number of time steps
  timestep_number = 0;

  // Counts total time  
  time = 0;
 
  // Here, we choose a time-stepping scheme that
  // is based on finite differences:
  // BE         = backward Euler scheme 
  // CN         = Crank-Nicolson scheme
  // CN_shifted = time-shifted Crank-Nicolson scheme 
  // For further properties of these schemes,
  // we refer to standard literature.
  if (time_stepping_scheme == "BE")
    theta = 1.0;
  else if (time_stepping_scheme == "CN")
    theta = 0.5;
  else if (time_stepping_scheme == "CN_shifted")
    theta = 0.5 + timestep;
  else 
    std::cout << "No such timestepping scheme" << std::endl;

  // In the following, we read a *.inp grid from a file.
  std::string grid_name;
  grid_name  = "Schraube_2d.msh"; 
 
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==2, ExcInternalError());
  grid_in.read_msh (input_file); 

  typename Triangulation<dim>::cell_iterator
    cell = triangulation.begin(),
    endc = triangulation.end();
    for (; cell!=endc; ++cell)
      for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
	if (cell->at_boundary(i))
	  if (cell->face(i)->center()(1) == 0)
	    {
	      cell->face(i)->set_boundary_id(3);
	    }
	  else if (cell->face(i)->center()(1) == -17.40)
	    {
	      cell->face(i)->set_boundary_id(2);
	    }
    
  
  global_refinement_steps = 1;
  pred_corr_levels = 0;
  triangulation.refine_global (global_refinement_steps); 

  filename_basis  = "solution_screw_domi_"; 
  //filename_basis  = "solution_test_"; 
  bool_use_error_oriented_Newton = false;
  bool_use_modified_Newton = true; // if true and be used, then must set error oriented false
  bool_set_explicitely_delta_fp = false; // if true, must set use_modified Newton to false
  bool_set_explicitely_constant_k = false;
  bool_use_pf_extra = false;

  a_fp = 0.01; // 0.01 for modified Newton
  b_fp = 5.0; // 5.0 for modified Newton
  max_no_of_augmented_L_penal_iterations = 50; // 50
  if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

 tolerance_augmented_L_iteration = 1.0e-3; //1.0e-4;
  tolerance_absolute_augmented_L_iteration = 1.0e-3;

  max_no_newton_steps  = 100;
  max_no_line_search_steps = 10;
 if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  lower_bound_newton_update = 1.0e-8;
  lower_bound_newton_residuum = 1.0e-8; 

  // When `true' be careful because tolerances are fixed 
  //  and differ from the choices here.
  bool_use_adaptive_newton_bound = false;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();

  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = false;

}


template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_Sneddon3D ()
{

  // Parameters
  current_pressure = 1.0e-3; 
  alpha_biot = 0.0;

  G_c = 1.0;
  delta_penal = 0.0; // simple penalization
  gamma_penal = 1.0e+2; //1.0e+1;  // augmented Lagrangian penalization

  gravity_x = 0.0;
  gravity_y = 0.0;
  volume_source = 0.0;


  // Elasticity parameters
  force_structure_x = 0.0;
  force_structure_y = 0.0;

  // Traction 
  traction_x = 0.0;
  traction_y = 0.0; 

  // Structure parameters
  density_structure = 1.0; 
  poisson_ratio_nu = 0.2; 
  E_modulus = 1.0;
  E_prime = E_modulus/(1.0 - poisson_ratio_nu * poisson_ratio_nu);
  
  lame_coefficient_mu = E_modulus/(2.0*(1 + poisson_ratio_nu));
  
  lame_coefficient_lambda = (2 * poisson_ratio_nu * lame_coefficient_mu)/(1.0 - 2 * poisson_ratio_nu);

  
  // Timestepping schemes
  //BE, CN, CN_shifted
  time_stepping_scheme = "BE";

  // Timestep size:
  timestep = 1.0; //2.5e-1; //1.0;

  // Maximum number of timesteps:
  max_no_timesteps = 5;
  end_time_value = 1e+10;
  
  // A variable to count the number of time steps
  timestep_number = 0;

  // Counts total time  
  time = 0;
 
  // Here, we choose a time-stepping scheme that
  // is based on finite differences:
  // BE         = backward Euler scheme 
  // CN         = Crank-Nicolson scheme
  // CN_shifted = time-shifted Crank-Nicolson scheme 
  // For further properties of these schemes,
  // we refer to standard literature.
  if (time_stepping_scheme == "BE")
    theta = 1.0;
  else if (time_stepping_scheme == "CN")
    theta = 0.5;
  else if (time_stepping_scheme == "CN_shifted")
    theta = 0.5 + timestep;
  else 
    std::cout << "No such timestepping scheme" << std::endl;

  // In the following, we read a *.inp grid from a file.
  std::string grid_name;
  grid_name  = "unit_cube_10.inp"; 
  
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==3, ExcInternalError());
  grid_in.read_ucd (input_file); 


 // TODO: insert all time step tolerances etc.
  
  global_refinement_steps = 2;
  pred_corr_levels = 0;
  triangulation.refine_global (global_refinement_steps); 

  //filename_basis  = "solution_pressurized_pf_extra_eps_2h_ref_7_adaptive_"; 
  filename_basis  = "solution_Sneddon_3D_"; 
  filename_basis_cod = "cod_ref_6_";
  bool_use_error_oriented_Newton = false;
  bool_use_pf_extra = false;
  max_no_of_augmented_L_penal_iterations = 10; //15; //40;
 if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

  tolerance_augmented_L_iteration = 1.0e-3; // 1e-4;
  tolerance_absolute_augmented_L_iteration = 1.0e-8;

  max_no_newton_steps  = 20; //20; // 50
  max_no_line_search_steps = 10; 
  if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  lower_bound_newton_update = 1.0e-8;
  lower_bound_newton_residuum = 1.0e-8; 

  bool_use_adaptive_newton_bound = false;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();

  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = false;

}



template <int dim>
void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_Het3D ()
{

  // Parameters
  current_pressure = 1.0e-3; 
  alpha_biot = 0.0;

  G_c = 1.0;
  delta_penal = 0.0; // simple penalization
  gamma_penal = 1.0e+2; //1.0e+1;  // augmented Lagrangian penalization

  gravity_x = 0.0;
  gravity_y = 0.0;
  volume_source = 0.0;


  // Elasticity parameters
  force_structure_x = 0.0;
  force_structure_y = 0.0;

  // Traction 
  traction_x = 0.0;
  traction_y = 0.0; 

  // Structure parameters
  density_structure = 1.0; 
  poisson_ratio_nu = 0.2; 
  E_modulus = 1.0;
  E_prime = E_modulus/(1.0 - poisson_ratio_nu * poisson_ratio_nu);
  
  lame_coefficient_mu = E_modulus/(2.0*(1 + poisson_ratio_nu));
  
  lame_coefficient_lambda = (2 * poisson_ratio_nu * lame_coefficient_mu)/(1.0 - 2 * poisson_ratio_nu);

  
  // Timestepping schemes
  //BE, CN, CN_shifted
  time_stepping_scheme = "BE";

  // Timestep size:
  timestep = 1.0; //2.5e-1; //1.0;

  // Maximum number of timesteps:
  max_no_timesteps = 30;
  end_time_value = 1e+10;
  
  // A variable to count the number of time steps
  timestep_number = 0;

  // Counts total time  
  time = 0;
 
  // Here, we choose a time-stepping scheme that
  // is based on finite differences:
  // BE         = backward Euler scheme 
  // CN         = Crank-Nicolson scheme
  // CN_shifted = time-shifted Crank-Nicolson scheme 
  // For further properties of these schemes,
  // we refer to standard literature.
  if (time_stepping_scheme == "BE")
    theta = 1.0;
  else if (time_stepping_scheme == "CN")
    theta = 0.5;
  else if (time_stepping_scheme == "CN_shifted")
    theta = 0.5 + timestep;
  else 
    std::cout << "No such timestepping scheme" << std::endl;

  // In the following, we read a *.inp grid from a file.
  std::string grid_name;
  grid_name  = "unit_cube_10.inp"; 
  
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file(grid_name.c_str());      
  Assert (dim==3, ExcInternalError());
  grid_in.read_ucd (input_file); 


 // TODO: insert all time step tolerances etc.
  
  global_refinement_steps = 2;
  pred_corr_levels = 0;
  triangulation.refine_global (global_refinement_steps);

  // TODO check E_modulus value in prm.file
  func_emodulus = new BitmapFunction<dim>("test.pgm",0,10,0,10,E_modulus,10.0*E_modulus);


  //filename_basis  = "solution_pressurized_pf_extra_eps_2h_ref_7_adaptive_"; 
  filename_basis  = "solution_Het_3D_a_001_5_"; 
  filename_basis_cod = "cod_ref_6_";
  bool_use_error_oriented_Newton = false;
  bool_use_modified_Newton = true; // if true and be used, then must set error oriented false
  bool_set_explicitely_delta_fp = false; // if true, must set use_modified Newton to false
  bool_set_explicitely_constant_k = false;
  bool_use_pf_extra = false;

  a_fp = 0.01; 
  b_fp = 5; 

  max_no_of_augmented_L_penal_iterations = 10; //15; //40;
  if (bool_use_strain_history)
    max_no_of_augmented_L_penal_iterations = 0;

  tolerance_augmented_L_iteration = 1.0e-3; // 1e-4;
  tolerance_absolute_augmented_L_iteration = 1.0e-8;

  max_no_newton_steps  = 100; //20; // 50
  max_no_line_search_steps = 10; 
  if (bool_use_modified_Newton)
    max_no_line_search_steps = 0;

  lower_bound_newton_update = 1.0e-8;
  lower_bound_newton_residuum = 1.0e-8; 

  bool_use_adaptive_newton_bound = true;
  // Required for heterogeneous materials
  //make random material vectors
  //make_material_vectors ();

  bool_plot_additional_solutions = false;
  bool_use_stress_splitting = false;

}







template <int dim>
void Dynamic_Fracture_Problem<dim>::make_material_vectors ()
{
  
  // tmp vector is used to write everything in one single vector
  // Later, this is split into 4 different vectors
  // We need 5 arguments because in material_ids are in the first column
  int no_arg = 3; // cell + mu + lame 
  Vector<double> tmp_vector(no_arg * triangulation.n_active_cells());

  // 2 vectors for the specific information on each cell
  lame_coefficient_mu_vector.reinit(triangulation.n_active_cells());
  lame_coefficient_lambda_vector.reinit(triangulation.n_active_cells());


  // Write everything from the single vector 
  // into the specific vector
  unsigned int l = 4;
  for (unsigned int i=0; i<triangulation.n_active_cells()/l; i++)
    {
      double factor = std::rand() % 500000; 
      for (unsigned int k=0; k<l; k++)
	{
	  lame_coefficient_mu_vector(l*i + k)   = factor + lame_coefficient_mu;
	  lame_coefficient_lambda_vector(l*i + k) = factor + lame_coefficient_lambda;
	}
    }
  
}
 




// This function is similar to many deal.II tuturial steps.
template <int dim>
void Dynamic_Fracture_Problem<dim>::setup_system ()
{
  timer.enter_section("Setup system.");


  value_phase_field_for_refinement = 0.9;


  // We set runtime parameters to drive the problem.
  // These parameters could also be read from a parameter file that
  // can be handled by the ParameterHandler object (see step-19)
  system_matrix.clear ();
  
  dof_handler.distribute_dofs (fe);  
  DoFRenumbering::Cuthill_McKee (dof_handler);

  std::vector<unsigned int> block_component (dim+1+dim,0);
  block_component[dim] = 1;
  block_component[dim+1] = 2;
  block_component[dim+2] = 2;
  DoFRenumbering::component_wise (dof_handler, block_component);

  {				 
    constraints.clear ();
    set_newton_bc ();
    DoFTools::make_hanging_node_constraints (dof_handler,
					     constraints);
  }
  constraints.close ();
  
  std::vector<unsigned int> dofs_per_block (3);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);  
  const unsigned int n_u = dofs_per_block[0],
    n_c =  dofs_per_block[1], n_v = dofs_per_block[2];

  std::cout << "Cells:\t"
            << triangulation.n_active_cells()
            << std::endl  	  
            << "DoFs:\t"
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_c << '+' << n_v <<  ')'
            << std::endl;


 
      
 {
   BlockDynamicSparsityPattern csp (3,3);

    csp.block(0,0).reinit (n_u, n_u);
    csp.block(0,1).reinit (n_u, n_c);
    csp.block(0,2).reinit (n_u, n_v);
  
    csp.block(1,0).reinit (n_c, n_u);
    csp.block(1,1).reinit (n_c, n_c);
    csp.block(1,2).reinit (n_c, n_v);

    csp.block(2,0).reinit (n_v, n_u);
    csp.block(2,1).reinit (n_v, n_c);
    csp.block(2,2).reinit (n_v, n_v);

 
    csp.collect_sizes();    
  

    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);

    sparsity_pattern.copy_from (csp);
  }
 
 system_matrix.reinit (sparsity_pattern);

  // Actual solution at time step n
  solution.reinit (3);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_c);
  solution.block(2).reinit (n_v);
  
  solution.collect_sizes ();
 
  // Old timestep solution at time step n-1
  old_timestep_solution.reinit (3);
  old_timestep_solution.block(0).reinit (n_u);
  old_timestep_solution.block(1).reinit (n_c);
  old_timestep_solution.block(2).reinit (n_v);
 
  old_timestep_solution.collect_sizes ();

  // Old timestep solution at time step n-2
  old_old_timestep_solution.reinit (3);
  old_old_timestep_solution.block(0).reinit (n_u);
  old_old_timestep_solution.block(1).reinit (n_c);
  old_old_timestep_solution.block(2).reinit (n_v);
 
  old_old_timestep_solution.collect_sizes ();

  // temporary solution for time adaptivity 
  tmp_solution_for_time_adaptivity.reinit (3);
  tmp_solution_for_time_adaptivity.block(0).reinit (n_u);
  tmp_solution_for_time_adaptivity.block(1).reinit (n_c);
  tmp_solution_for_time_adaptivity.block(2).reinit (n_v);

  tmp_solution_for_time_adaptivity.collect_sizes ();


  // Updates for Newton's method
  newton_update.reinit (3);
  newton_update.block(0).reinit (n_u);
  newton_update.block(1).reinit (n_c);
  newton_update.block(2).reinit (n_v);

  newton_update.collect_sizes ();
 
  // Residual for  Newton's method
  system_rhs.reinit (3);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_c);
  system_rhs.block(2).reinit (n_v);

  system_rhs.collect_sizes ();


  // Lambda penal function
  solution_lambda_penal_func.reinit (3);
  solution_lambda_penal_func.block(0).reinit (n_u);
  solution_lambda_penal_func.block(1).reinit (n_c);
  solution_lambda_penal_func.block(2).reinit (n_v);

  solution_lambda_penal_func.collect_sizes ();

  old_timestep_solution_lambda_penal_func.reinit (3);
  old_timestep_solution_lambda_penal_func.block(0).reinit (n_u);
  old_timestep_solution_lambda_penal_func.block(1).reinit (n_c);
  old_timestep_solution_lambda_penal_func.block(2).reinit (n_v);

  old_timestep_solution_lambda_penal_func.collect_sizes ();


  // TODO: check, if this needs to be updated
  // after refining the mesh (currently it is because
  // setup_system is called in the refine_mesh function)
  if (bool_use_strain_history)
    setup_quadrature_point_history ();

  double min_cell_diameter = 1.0e+10;
 typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  for (; cell!=endc; ++cell)
    { 
      cell_diameter = cell->diameter();
      if (min_cell_diameter > cell_diameter)
	min_cell_diameter = cell_diameter;	

    }

  std::cout << "Min cell dia: " << min_cell_diameter << std::endl;

  timer.exit_section(); 
}





// In this function, we assemble the Jacobian matrix
// for the Newton iteration. 
template <int dim>
void Dynamic_Fracture_Problem<dim>::assemble_system_matrix ()
{
  timer.enter_section("Assemble Matrix.");
  system_matrix=0;
     
  QGauss<dim>   quadrature_formula(degree+2);  
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values |
                           update_gradients);
  
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
				    update_values         | update_quadrature_points  |
				    update_normal_vectors | update_gradients |
				    update_JxW_values);
   
  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  
  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int n_face_q_points   = face_quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell); 
		

  // Now, we are going to use the 
  // FEValuesExtractors to determine
  // the four principle variables
  const FEValuesExtractors::Vector displacements (0); // 2
  const FEValuesExtractors::Scalar phase_field (dim); // 4
  const FEValuesExtractors::Vector velocities (dim+1); // 4

  // We declare Vectors and Tensors for 
  // the solutions at the previous Newton iteration:
  std::vector<Vector<double> > old_solution_values (n_q_points, 
				 		    Vector<double>(dim+1+dim));

  std::vector<std::vector<Tensor<1,dim> > > old_solution_grads (n_q_points, 
								std::vector<Tensor<1,dim> > (dim+1+dim));

  std::vector<Vector<double> >  old_solution_face_values (n_face_q_points, 
							  Vector<double>(dim+1+dim));
       
  std::vector<std::vector<Tensor<1,dim> > > old_solution_face_grads (n_face_q_points, 
								     std::vector<Tensor<1,dim> > (dim+1+dim));
    

  // We declare Vectors and Tensors for 
  // the solution at the previous time step:
   std::vector<Vector<double> > old_timestep_solution_values (n_q_points, 
				 		    Vector<double>(dim+1+dim));

   std::vector<Vector<double> >   old_old_timestep_solution_values (n_q_points, 
								    Vector<double>(dim+1+dim));


  std::vector<std::vector<Tensor<1,dim> > > old_timestep_solution_grads (n_q_points, 
  					  std::vector<Tensor<1,dim> > (dim+1+dim));


  std::vector<Vector<double> >   old_timestep_solution_face_values (n_face_q_points, 
								    Vector<double>(dim+1+dim));
  
    
  std::vector<std::vector<Tensor<1,dim> > >  old_timestep_solution_face_grads (n_face_q_points, 
									       std::vector<Tensor<1,dim> > (dim+1+dim));

std::vector<Vector<double> > old_solution_values_lambda_penal_func (n_q_points, 
								     Vector<double>(dim+1+dim));

   
  // Declaring test functions:
  std::vector<Tensor<1,dim> > phi_i_u (dofs_per_cell); 
  std::vector<Tensor<2,dim> > phi_i_grads_u(dofs_per_cell);
  std::vector<double>         phi_i_pf(dofs_per_cell); 
  std::vector<Tensor<1,dim> > phi_i_grads_pf (dofs_per_cell);
  // u and v here components of displacements
  std::vector<Tensor<1,dim> > phi_i_v (dofs_per_cell); 
  std::vector<Tensor<2,dim> > phi_i_grads_v(dofs_per_cell);

  // This is the identity matrix in two dimensions:
  const Tensor<2,dim> Identity = ALE_Transformations
    ::get_Identity<dim> ();

  Tensor<2,dim> zero_matrix;
  zero_matrix.clear();
 				     				   
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  unsigned int cell_counter = 0;
  for (; cell!=endc; ++cell)
    { 
      //lame_coefficient_mu = lame_coefficient_mu_vector(cell_counter);
      //lame_coefficient_lambda = lame_coefficient_lambda_vector(cell_counter);
      cell_counter ++;

      fe_values.reinit (cell);
      local_matrix = 0;
      
      // We need the cell diameter to control the fluid mesh motion
      cell_diameter = cell->diameter();

      if (test_case == "Het3D")
	{
	  E_modulus = func_emodulus->value(cell->center(), 0);
	  E_modulus += 1.0;
	  //E_modulus = 5.0;
	  
	  lame_coefficient_mu = E_modulus / (2.0 * (1 + poisson_ratio_nu));
	  
	  lame_coefficient_lambda = (2 * poisson_ratio_nu * lame_coefficient_mu)
	    / (1.0 - 2 * poisson_ratio_nu);
	}

      // Old Newton iteration values
      fe_values.get_function_values (solution, old_solution_values);
      fe_values.get_function_gradients (solution, old_solution_grads);
      
      // Old_timestep_solution values
      fe_values.get_function_values (old_timestep_solution, old_timestep_solution_values);
      fe_values.get_function_gradients (old_timestep_solution, old_timestep_solution_grads);

      // Old_Old_timestep_solution values
      fe_values.get_function_values (old_old_timestep_solution, old_old_timestep_solution_values);
  
   // Old Newton iteration values lambda 
      fe_values.get_function_values (solution_lambda_penal_func, old_solution_values_lambda_penal_func);
   
  
  const PointHistory<dim> *local_quadrature_points_data
	= reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

	  for (unsigned int q=0; q<n_q_points; ++q)
	    {	      
	      for (unsigned int k=0; k<dofs_per_cell; ++k)
		{
		  phi_i_u[k]       = fe_values[displacements].value (k, q);
		  phi_i_grads_u[k] = fe_values[displacements].gradient (k, q);
		  phi_i_pf[k]       = fe_values[phase_field].value (k, q);
		  phi_i_grads_pf[k] = fe_values[phase_field].gradient (k, q);
		  phi_i_v[k]       = fe_values[velocities].value (k, q);
		  phi_i_grads_v[k] = fe_values[velocities].gradient (k, q);
		}
	      

	      const double pf = old_solution_values[q](dim); 
	      const double old_timestep_pf = old_timestep_solution_values[q](dim); 
	      const double old_old_timestep_pf = old_old_timestep_solution_values[q](dim); 

	      const double  lambda_penal_func = old_solution_values_lambda_penal_func[q](dim);

	      double pf_extra = pf;
              // Linearization by extrapolation to cope with non-convexity of the underlying
              // energy functional.
              // This idea might be refined in a future work (be also careful because
              // theoretically, we do not have time regularity; therefore extrapolation in time
              // might be questionable. But for the time being, this is numerically robust.
              pf_extra = old_old_timestep_pf + (time - (time-old_timestep-old_old_timestep))/
                         (time-old_timestep - (time-old_timestep-old_old_timestep)) * (old_timestep_pf - old_old_timestep_pf);


              if (pf_extra <= 0.0)
                pf_extra = 0.0;
              if (pf_extra >= 1.0)
                pf_extra = 1.0;

	      const Tensor<2,dim> grad_u = ALE_Transformations 
		::get_grad_u<dim> (q, old_solution_grads);

	      double divergence_u = old_solution_grads[q][0][0] +  old_solution_grads[q][1][1];
	      if (dim == 3)
		divergence_u += old_solution_grads[q][2][2];

	      // Linearized strain
	      const Tensor<2,dim> E = 0.5 * (grad_u + transpose(grad_u));
	      
	      const double tr_E = Structure_Terms_in_ALE
		::get_tr_E<dim> (E);

	      
	      Tensor<2,dim> stress_term;
	      stress_term.clear();
	      stress_term = lame_coefficient_lambda * tr_E * Identity
		+ 2 * lame_coefficient_mu * E;


	      Tensor<2,dim> stress_term_plus;
              Tensor<2,dim> stress_term_minus;
	      
	      if ((timestep_number > 1) && (dim == 2) && bool_use_stress_splitting)
		{
		  decompose_stress(stress_term_plus, stress_term_minus,
                                   E, tr_E, zero_matrix , 0.0,
                                   lame_coefficient_lambda,
                                   lame_coefficient_mu, false);
		}
	      else 
		{
		  stress_term_plus = stress_term;
		  stress_term_minus = 0;
		}

	      const Tensor<2,dim> &stress_term_history
		= local_quadrature_points_data[q].old_stress;

	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  // Simple penalization
//		  double pf_minus_old_timestep_pf_plus = 0.0;
//		  if ((pf - old_timestep_pf) < 0.0)
//		    pf_minus_old_timestep_pf_plus = 0.0;
//		  else 
//		    pf_minus_old_timestep_pf_plus = phi_i_pf[i]; 


		  double chi = 0.0;
		  if ((lambda_penal_func + gamma_penal * (pf - old_timestep_pf)) > 0.0)
		    chi = 1.0;
		  else 
		    chi = 0.0;

		  double divergence_u_LinU = phi_i_grads_u[i][0][0] + phi_i_grads_u[i][1][1];
		  if (dim == 3)
		    divergence_u_LinU += phi_i_grads_u[i][2][2];
	    	     
		
		  const Tensor<2, dim> E_LinU = 0.5
		    * (phi_i_grads_u[i] + transpose(phi_i_grads_u[i]));

		  const double tr_E_LinU = dealii::trace(E_LinU);

		  Tensor<2,dim> stress_term_LinU;
		  stress_term_LinU = lame_coefficient_lambda * tr_E_LinU * Identity
		    + 2 * lame_coefficient_mu * E_LinU;

		  Tensor<2,dim> stress_term_plus_LinU;
		  Tensor<2,dim> stress_term_minus_LinU;
		  if ((timestep_number > 1) && (dim == 2) && bool_use_stress_splitting)
		    {
		  decompose_stress(stress_term_plus_LinU, stress_term_minus_LinU,
				   E, tr_E, E_LinU, tr_E_LinU,
				   lame_coefficient_lambda,
				   lame_coefficient_mu,
				   true);
		    }
		  else 
		    {
		      stress_term_plus_LinU = stress_term_LinU;
		      stress_term_minus_LinU = 0;
		    }


		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    {
		      // STVK 
		      const unsigned int comp_j = fe.system_to_component_index(j).first; 
		      if (comp_j < dim)
			{
			  if (bool_use_dynamic_code_with_velocities)
			    {
			      // TODO: density, dividing over time step
			      local_matrix(j,i) += phi_i_v[i] * phi_i_u[j] * fe_values.JxW(q);
			    }

			  // pf is solution variable
			  if (!bool_use_pf_extra)
			    {
			      // Simple Hessian modification (not very reliable)
			      //local_matrix(j,i) += 0.0 * timestep * theta * 
			      //	1.0/(cell_diameter*cell_diameter) * 0.1 * lame_coefficient_mu * (phi_i_u[i] * phi_i_u[j]
			      //	 ) * fe_values.JxW(q);



			      local_matrix(j,i) += timestep * theta * 
				(delta_fixed_point_newton * dealii::scalar_product((1-constant_k) * 2.0 * pf * phi_i_pf[i] * stress_term_plus,phi_i_grads_u[j])
				 + dealii::scalar_product(((1-constant_k) * pf * pf + constant_k) * stress_term_plus_LinU, phi_i_grads_u[j])
				 + dealii::scalar_product(stress_term_minus_LinU, phi_i_grads_u[j])
				 // Pressure (pf is solution variable)
				 - 2.0 * (alpha_biot - 1.0) * current_pressure * pf * phi_i_pf[i] * divergence_u_LinU
				 ) * fe_values.JxW(q);
			    }
			  else if (bool_use_pf_extra)
			    {
			      // pf extrapolated
			      local_matrix(j,i) += timestep * theta * 
				(dealii::scalar_product(((1-constant_k) * pf_extra * pf_extra + constant_k) * 	      
						stress_term_plus_LinU, phi_i_grads_u[j])
				 + dealii::scalar_product(stress_term_minus_LinU, phi_i_grads_u[j])
				 ) * fe_values.JxW(q);
			    }
			}		     
		      else if (comp_j == dim)
			{

			  if (!bool_use_strain_history)
			    {

			  // Simple penalization
			  //local_matrix(j,i) += delta_penal *  1.0/cell_diameter *  phi_i_pf[j] 
			  //  * pf_minus_old_timestep_pf_plus* fe_values.JxW(q);

			  // Augmented Lagrangian penalization
			      local_matrix(j,i) += chi * gamma_penal * phi_i_pf[i] * phi_i_pf[j] * fe_values.JxW(q);


			      local_matrix(j,i) += timestep * theta *
				( (1-constant_k) * (dealii::scalar_product(stress_term_plus_LinU, E)
						    + dealii::scalar_product(stress_term_plus, E_LinU)) * pf * phi_i_pf[j]
				  +(1-constant_k) * dealii::scalar_product(stress_term_plus, E) * phi_i_pf[i] * phi_i_pf[j]
				  + G_c/alpha_eps * phi_i_pf[i] * phi_i_pf[j]  
				  + G_c * alpha_eps * phi_i_grads_pf[i] * phi_i_grads_pf[j]
				  // Pressure terms
				  - 2.0 * (alpha_biot - 1.0) * current_pressure *
				  (pf * divergence_u_LinU + phi_i_pf[i] * divergence_u) * phi_i_pf[j]
				  ) * fe_values.JxW(q);      

			    }
			  else if (bool_use_strain_history)
			    {

			      local_matrix(j,i) += timestep * theta *
			    ((1-constant_k) * (0.0 * dealii::scalar_product(stress_term_plus_LinU, E)
			     		       + dealii::scalar_product(stress_term_history, E_LinU)) * pf * phi_i_pf[j]
			     +(1-constant_k) * dealii::scalar_product(stress_term_history, E) * phi_i_pf[i] * phi_i_pf[j]

			     + G_c/alpha_eps * phi_i_pf[i] * phi_i_pf[j]  
			     + G_c * alpha_eps * phi_i_grads_pf[i] * phi_i_grads_pf[j]
			     // Pressure terms
			     - 2.0 * (alpha_biot - 1.0) * current_pressure *
			     (pf * divergence_u_LinU + phi_i_pf[i] * divergence_u) * phi_i_pf[j]
			     ) * fe_values.JxW(q);     

			    }

			}
		      else if (comp_j > dim)
			{
			  //local_matrix(j,i) += scalar_product(phi_i_grads_v[i],phi_i_grads_v[j])
			  // * fe_values.JxW(q);

			  if (bool_use_dynamic_code_with_velocities)
			    {
			      local_matrix(j,i) += (phi_i_u[i] - timestep * theta * phi_i_v[i])  * phi_i_v[j] * fe_values.JxW(q);
			    }
			  else 
			    {
			      local_matrix(j,i) += phi_i_v[i]  * phi_i_v[j]
				* fe_values.JxW(q);
			    }


			}

		      
		    }  // end j dofs
		  	     
		}  // end i dofs 
	      
	    }   // end n_q_points  

	  
	  cell->get_dof_indices (local_dof_indices);
	  constraints.distribute_local_to_global (local_matrix, local_dof_indices,
						  system_matrix);

      // end cell
    }   
  
  timer.exit_section();
}



// In this function we assemble the semi-linear 
// of the right hand side of Newton's method (its residual).
// The framework is in principal the same as for the 
// system matrix.
template <int dim>
void
Dynamic_Fracture_Problem<dim>::assemble_system_rhs ()
{
  timer.enter_section("Assemble Rhs.");
  system_rhs=0;
  
  QGauss<dim>   quadrature_formula(degree+2);
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values |
                           update_gradients);

  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
				    update_values         | update_quadrature_points  |
				    update_normal_vectors | update_gradients |
				    update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  
  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int n_face_q_points   = face_quadrature_formula.size();
 
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
  const FEValuesExtractors::Vector displacements (0); 
  const FEValuesExtractors::Scalar phase_field (dim);
  const FEValuesExtractors::Vector velocities (dim+1);
 
  std::vector<Vector<double> > 
    old_solution_values (n_q_points, Vector<double>(dim+1+dim));

  std::vector<std::vector<Tensor<1,dim> > > 
    old_solution_grads (n_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));


  std::vector<Vector<double> > 
    old_solution_face_values (n_face_q_points, Vector<double>(dim+1+dim));
  
  std::vector<std::vector<Tensor<1,dim> > > 
    old_solution_face_grads (n_face_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));
  
  std::vector<Vector<double> > 
    old_timestep_solution_values (n_q_points, Vector<double>(dim+1+dim));

  std::vector<Vector<double> > 
    old_old_timestep_solution_values (n_q_points, Vector<double>(dim+1+dim));

  std::vector<std::vector<Tensor<1,dim> > > 
    old_timestep_solution_grads (n_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));

  std::vector<Vector<double> > 
    old_timestep_solution_face_values (n_face_q_points, Vector<double>(dim+1+dim));
     
  std::vector<std::vector<Tensor<1,dim> > > 
    old_timestep_solution_face_grads (n_face_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));

 std::vector<Vector<double> > old_solution_values_lambda_penal_func (n_q_points, 
								     Vector<double>(dim+1+dim));
 
   
  const Tensor<2,dim> Identity = ALE_Transformations
    ::get_Identity<dim> ();

  Tensor<2,dim> zero_matrix;
  zero_matrix.clear();

  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  unsigned int cell_counter = 0;
  for (; cell!=endc; ++cell)
    { 
      //Required for heterogeneous media
      //lame_coefficient_mu = lame_coefficient_mu_vector(cell_counter);
      //lame_coefficient_lambda = lame_coefficient_lambda_vector(cell_counter);
      cell_counter++;


      fe_values.reinit (cell);	 
      local_rhs = 0;   	
      
      cell_diameter = cell->diameter();

      if (test_case == "Het3D")
	{
	  E_modulus = func_emodulus->value(cell->center(), 0);
	  E_modulus += 1.0;
	  //E_modulus = 5.0;
	  
	  lame_coefficient_mu = E_modulus / (2.0 * (1 + poisson_ratio_nu));
	  
	  lame_coefficient_lambda = (2 * poisson_ratio_nu * lame_coefficient_mu)
	    / (1.0 - 2 * poisson_ratio_nu);
	}

      
      // old Newton iteration
      fe_values.get_function_values (solution, old_solution_values);
      fe_values.get_function_gradients (solution, old_solution_grads);
            
      // old timestep iteration
      fe_values.get_function_values (old_timestep_solution, old_timestep_solution_values);
      fe_values.get_function_gradients (old_timestep_solution, old_timestep_solution_grads);

      // old timestep iteration
      fe_values.get_function_values (old_old_timestep_solution, old_old_timestep_solution_values);

      // Old Newton iteration values lambda 
      fe_values.get_function_values (solution_lambda_penal_func, old_solution_values_lambda_penal_func);



    const PointHistory<dim> *local_quadrature_points_data
	= reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

      

	  for (unsigned int q=0; q<n_q_points; ++q)
	    {
	      const double pf = old_solution_values[q](dim);
	      const double old_timestep_pf = old_timestep_solution_values[q](dim); 

	      const double old_old_timestep_pf = old_old_timestep_solution_values[q](dim); 
	      
	      const double  lambda_penal_func = old_solution_values_lambda_penal_func[q](dim);

	      double pf_extra = pf;
              // Linearization by extrapolation to cope with non-convexity of the underlying
              // energy functional.
              // This idea might be refined in a future work (be also careful because
              // theoretically, we do not have time regularity; therefore extrapolation in time
              // might be questionable. But for the time being, this is numerically robust.
              pf_extra = old_old_timestep_pf + (time - (time-old_timestep-old_old_timestep))/
                         (time-old_timestep - (time-old_timestep-old_old_timestep)) * (old_timestep_pf - old_old_timestep_pf);


              if (pf_extra <= 0.0)
                pf_extra = 0.0;
              if (pf_extra >= 1.0)
                pf_extra = 1.0;


	      // Simple penalization
//	      double pf_minus_old_timestep_pf_plus = 0.0;
//	      if ((pf - old_timestep_pf) < 0.0)
//		pf_minus_old_timestep_pf_plus = 0.0;
//	      else 
//		pf_minus_old_timestep_pf_plus = pf - old_timestep_pf;


	      // Augmented Lagrangian; but first term missing
	      double pf_minus_old_timestep_pf_plus = 0.0;
	      if ((lambda_penal_func + gamma_penal * (pf - old_timestep_pf)) < 0.0)
		pf_minus_old_timestep_pf_plus = 0.0;
	      else 
		pf_minus_old_timestep_pf_plus = lambda_penal_func + gamma_penal * (pf - old_timestep_pf);

	      //std::cout << pf_minus_old_timestep_pf_plus << std::endl;

	      const Tensor<1,dim> grad_pf = ALE_Transformations
		::get_grad_pf<dim> (q, old_solution_grads);

	      const Tensor<1,dim> old_timestep_grad_pf = ALE_Transformations
		::get_grad_pf<dim> (q, old_timestep_solution_grads);

	      const Tensor<1,dim> v = ALE_Transformations 
		::get_v<dim> (q, old_solution_values);

	      const Tensor<1,dim> old_timestep_v = ALE_Transformations 
		::get_v<dim> (q, old_timestep_solution_values);

	      const Tensor<2,dim> grad_v = ALE_Transformations 
		::get_grad_v<dim> (q, old_solution_grads);
	      

	      const Tensor<1,dim> u = ALE_Transformations 
		::get_u<dim> (q, old_solution_values);
	      
	      const Tensor<1,dim> old_timestep_u = ALE_Transformations 
		::get_u<dim> (q, old_timestep_solution_values);
	      
	      const Tensor<2,dim> grad_u = ALE_Transformations 
		::get_grad_u<dim> (q, old_solution_grads);

	      const Tensor<2,dim> old_timestep_grad_u = ALE_Transformations 
		::get_grad_u<dim> (q, old_timestep_solution_grads);
	      
	      double divergence_u = old_solution_grads[q][0][0] +  old_solution_grads[q][1][1];
	      if (dim == 3)
		divergence_u += old_solution_grads[q][2][2];

	      double old_timestep_divergence_u = old_timestep_solution_grads[q][0][0] +  old_timestep_solution_grads[q][1][1];
	      if (dim == 3)
		old_timestep_divergence_u += old_timestep_solution_grads[q][2][2];

	      // Linearized strain
	      const Tensor<2,dim> E = 0.5 * (grad_u + transpose(grad_u));
	      
	      const double tr_E = Structure_Terms_in_ALE
		::get_tr_E<dim> (E);

	      const Tensor<2,dim> old_timestep_E = 0.5 * (old_timestep_grad_u + transpose(old_timestep_grad_u));
	      
	      const double old_timestep_tr_E = Structure_Terms_in_ALE
		::get_tr_E<dim> (old_timestep_E);
	      
	      
	      Tensor<2,dim> stress_term;
	      stress_term.clear();
	      stress_term = lame_coefficient_lambda * tr_E * Identity
		+ 2 * lame_coefficient_mu * E;

	      Tensor<2,dim> old_timestep_stress_term;
	      old_timestep_stress_term.clear();
	      old_timestep_stress_term = lame_coefficient_lambda * old_timestep_tr_E * Identity
		+ 2 * lame_coefficient_mu * old_timestep_E;

	      Tensor<2,dim> stress_term_plus;
              Tensor<2,dim> stress_term_minus;

	      Tensor<2,dim> old_timestep_stress_term_plus;
              Tensor<2,dim> old_timestep_stress_term_minus;
	      
	      if ((timestep_number > 1) && (dim == 2) && bool_use_stress_splitting)
		{
		  decompose_stress(stress_term_plus, stress_term_minus,
                                   E, tr_E, zero_matrix , 0.0,
                                   lame_coefficient_lambda,
                                   lame_coefficient_mu, false);

		  // Old timestep
		  decompose_stress(old_timestep_stress_term_plus, old_timestep_stress_term_minus,
                                   old_timestep_E, old_timestep_tr_E, zero_matrix , 0.0,
                                   lame_coefficient_lambda,
                                   lame_coefficient_mu, false);
		}
	      else 
		{
		  stress_term_plus = stress_term;
		  stress_term_minus = 0;

		  // Old timestep
		  old_timestep_stress_term_plus = old_timestep_stress_term;
		  old_timestep_stress_term_minus = 0;
		}
	      

	       const Tensor<2,dim> &stress_term_history
		= local_quadrature_points_data[q].old_stress; 


	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  // Info: Compute velocities
		  const unsigned int comp_i = fe.system_to_component_index(i).first; 
		  if (comp_i < dim)
		    { 
		      const Tensor<1,dim> phi_i_u = fe_values[displacements].value (i, q);
		      const Tensor<2,dim> phi_i_grads_u = fe_values[displacements].gradient (i, q);

		      double divergence_u_LinU = phi_i_grads_u[0][0] + phi_i_grads_u[1][1];
		      if (dim == 3)
			divergence_u_LinU += phi_i_grads_u[2][2];


		      // Acceleration term
		      if (bool_use_dynamic_code_with_velocities)
			{
			  local_rhs(i) -= (v - old_timestep_v) * phi_i_u * fe_values.JxW(q);
			}

		      // pf is solution variable
		      if (!bool_use_pf_extra)
			{
			  
			  // Current timestep solution
			  local_rhs(i) -= timestep * theta * 
			    (dealii::scalar_product(((1-constant_k) * pf * pf + constant_k) *	  
					    stress_term_plus, phi_i_grads_u)
			     + dealii::scalar_product(stress_term_minus, phi_i_grads_u)
			     // Pressure terms (pf is solution variable)
			     - (alpha_biot - 1.0) * current_pressure * pf * pf * divergence_u_LinU
			     ) * fe_values.JxW(q); 

			  // Old timestep solution
			  local_rhs(i) -= timestep * (1.0-theta) * 
			    (dealii::scalar_product(((1-constant_k) * old_timestep_pf * old_timestep_pf + constant_k) *	  
					    old_timestep_stress_term_plus, phi_i_grads_u)
			     + dealii::scalar_product(old_timestep_stress_term_minus, phi_i_grads_u)
			     // Pressure terms (pf is solution variable)
			     - (alpha_biot - 1.0) * old_timestep_current_pressure * old_timestep_pf * old_timestep_pf * divergence_u_LinU
			     ) * fe_values.JxW(q); 

			  



			}
		      else if (bool_use_pf_extra)
			{
			 
			  // Current timestep solution
			  // pf extrapolated
			  local_rhs(i) -= timestep * theta * 
			    (dealii::scalar_product(((1-constant_k) * pf_extra * pf_extra + constant_k) *	  
					    stress_term_plus, phi_i_grads_u)
			     + dealii::scalar_product(stress_term_minus, phi_i_grads_u)
			     // Pressure terms (extrapolated)
			     - (alpha_biot - 1.0) * current_pressure * pf_extra * pf_extra * divergence_u_LinU
			     ) * fe_values.JxW(q); 

			  // Old timestep solution
			  // TODO: need to define previous timestep extrapolation
			  local_rhs(i) -= timestep * (1.0 - theta) * 
			    (dealii::scalar_product(((1-constant_k) * pf_extra * pf_extra + constant_k) *	  
					    old_timestep_stress_term_plus, phi_i_grads_u)
			     + dealii::scalar_product(old_timestep_stress_term_minus, phi_i_grads_u)
			     // Pressure terms (extrapolated)
			     - (alpha_biot - 1.0) * old_timestep_current_pressure * pf_extra * pf_extra * divergence_u_LinU
			     ) * fe_values.JxW(q); 
			  
			  
			} 
		      
		    }		
		  else if (comp_i == dim)
		    {
		      const double phi_i_pf = fe_values[phase_field].value (i, q);
		      const Tensor<1,dim> phi_i_grads_pf = fe_values[phase_field].gradient (i, q);

		      if (!bool_use_strain_history)
			{
			  //  Simple penalization
			  //local_rhs(i) -= delta_penal *  1.0/cell_diameter  
			  //	* pf_minus_old_timestep_pf_plus * phi_i_pf * fe_values.JxW(q);

		     
			  // Augmented Lagrangian penalization
			  local_rhs(i) -= pf_minus_old_timestep_pf_plus * phi_i_pf * fe_values.JxW(q);
		      
			  // Current time step
			  local_rhs(i) -= timestep * theta * 
			    ((1.0 - constant_k) * dealii::scalar_product(stress_term_plus, E) * pf * phi_i_pf
			     - G_c/alpha_eps * (1.0 - pf) * phi_i_pf
			     + G_c * alpha_eps * grad_pf * phi_i_grads_pf
			     // Pressure terms
			     - 2.0 * (alpha_biot - 1.0) * current_pressure * pf * divergence_u * phi_i_pf
			     ) * fe_values.JxW(q);

			  // Old timestep
			  local_rhs(i) -= timestep * (1.0-theta) * 
			    ((1.0 - constant_k) * dealii::scalar_product(old_timestep_stress_term_plus, E) * old_timestep_pf * phi_i_pf
			     - G_c/alpha_eps * (1.0 - old_timestep_pf) * phi_i_pf
			     + G_c * alpha_eps * old_timestep_grad_pf * phi_i_grads_pf
			     // Pressure terms
			     - 2.0 * (alpha_biot - 1.0) * old_timestep_current_pressure * old_timestep_pf * old_timestep_divergence_u * phi_i_pf
			     ) * fe_values.JxW(q);
			  
			}
		      else if (bool_use_strain_history)
			{
			  std::cout << "Aborting ..." << std::endl;
			  abort();
			  // only implemented for BE (backward Euler) so far
//                         local_rhs(i) -= timestep * theta * 
//			((1.0 - constant_k) * scalar_product(stress_term_history, E) * pf * phi_i_pf
//			 - G_c/alpha_eps * (1.0 - pf) * phi_i_pf
//			 + G_c * alpha_eps * grad_pf * phi_i_grads_pf
//			 // Pressure terms
//			 - 2.0 * (alpha_biot - 1.0) * current_pressure * pf * divergence_u * phi_i_pf
//			 ) * fe_values.JxW(q);

			}

		    } 
		  else if (comp_i > dim)
		    { 
		      // Info: Compute displacements
		      const Tensor<1,dim> phi_i_v       = fe_values[velocities].value (i, q);
		      const Tensor<2,dim> phi_i_grads_v = fe_values[velocities].gradient (i, q);
		      
		      //local_rhs(i) -= 
		      //	(scalar_product(grad_v,phi_i_grads_v) - 1.0 * phi_i_v[0]
		      //	 ) * fe_values.JxW(q); 

		      if (bool_use_dynamic_code_with_velocities)
			{
			  local_rhs(i) -= ((u - old_timestep_u) - timestep * (theta * v + (1-theta) * old_timestep_v)) * phi_i_v * fe_values.JxW(q); 
			}
		      else 
			{
			  local_rhs(i) -= 
			    (v * phi_i_v
			     ) * fe_values.JxW(q); 
			}

		    }  

		  	  
		} // end i
	      		   
	    } // end n_q_points 
	  


	  // upper traction (non-pay zone)
	  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	    {
	      if (cell->face(face)->at_boundary() && 		  
		  (cell->face(face)->boundary_id() == 3) 
		  )
		{
		  
		  fe_face_values.reinit (cell, face);
		  
		  for (unsigned int q=0; q<n_face_q_points; ++q)
		    {	
		      Tensor<1,dim> neumann_value;
		      neumann_value[0] = 0.0; //time * traction_x;
		      neumann_value[1] = 0.0; //time * traction_y;
			
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
			  const unsigned int comp_i = fe.system_to_component_index(i).first; 
			  if (comp_i == 0 || comp_i == 1)
			    {  
			      local_rhs(i) +=  1.0 * (timestep * theta * 
						      neumann_value * fe_face_values[displacements].value (i, q) 
						      ) * fe_face_values.JxW(q);					   
			    }
			  // end i
			}  
		      // end face_n_q_points    
		    }                                     
		} 
	    }  // end face integrals upper fraction



	  // lower traction (non-pay zone)
	  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	    {
	      if (cell->face(face)->at_boundary() && 		  
		  (cell->face(face)->boundary_id() == 2) 
		  )
		{
		  
		  fe_face_values.reinit (cell, face);
		  
		  for (unsigned int q=0; q<n_face_q_points; ++q)
		    {	
		      Tensor<1,dim> neumann_value;
		      neumann_value[0] = - 0.0; //time * traction_x;
		      neumann_value[1] = - 0.0; //time * traction_y;
			
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
			  const unsigned int comp_i = fe.system_to_component_index(i).first; 
			  if (comp_i == 0 || comp_i == 1)
			    {  
			      local_rhs(i) +=  1.0 * (timestep * theta * 
						      neumann_value * fe_face_values[displacements].value (i, q) 
						      ) * fe_face_values.JxW(q);					   
			    }
			  // end i
			}  
		      // end face_n_q_points    
		    }                                     
		} 
	    }  // end face integrals lower traction










	  cell->get_dof_indices (local_dof_indices);
	  constraints.distribute_local_to_global (local_rhs, local_dof_indices,
						  system_rhs);
	  

      
    }  // end cell
      
  timer.exit_section();
}


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
	component_mask[0]     = false;
	component_mask[1]     = true;

	// changed component_mask[dim+1]     = false;
  component_mask[dim+1]     = false;
	component_mask[dim+2]     = true;
  // 0 is the left edge
	VectorTools::interpolate_boundary_values (dof_handler,
						  0,
						  ZeroFunction<dim>(dim+1+dim),
						  boundary_values,
						  component_mask);

	
  // 1 is the right edge
  component_mask[0]     = false;
	component_mask[1]     = true;
  // changed component_mask[dim+1]     = false;
	component_mask[dim+1]     = false;
	component_mask[dim+2]     = true;
  
	VectorTools::interpolate_boundary_values (dof_handler,
						  1,
						  ZeroFunction<dim>(dim+1+dim),
						  boundary_values,
						  component_mask);

	
  
  
  
  // 2 is the bottom edge
  component_mask[0]     = true;
	component_mask[1]     = true;
	component_mask[dim+1]     = true;
	component_mask[dim+2]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  2,
						  ZeroFunction<dim>(dim+1+dim),
						  boundary_values,
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
						  NonhomDirichletBoundaryValues<dim>(time, test_case,alpha_eps),
						  boundary_values,
						  component_mask);

	
  // 4 is the bottom part of the crack
  // uy=0, vy=0
  component_mask[0]     = false;
	component_mask[1]     = true;
	component_mask[dim+1]     = false;
	component_mask[dim+2]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  4,
						  ZeroFunction<dim>(dim+1+dim),
						  boundary_values,
						  component_mask);

	if (bool_initial_crack_via_phase_field)
	{
	// Phase-field
	component_mask[0]     = false;
	component_mask[1]     = false;
	component_mask[2]     = true;
	component_mask[dim+1]     = false;
	component_mask[dim+2]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  1,
						  NonhomDirichletBoundaryValues<dim>(time, test_case,alpha_eps),
						  boundary_values,
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
						  boundary_values,
						  component_mask);


	component_mask[0]     = false;
	component_mask[1]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  11,
						  ZeroFunction<dim>(dim+1+dim),
						  boundary_values,
						  component_mask);


	component_mask[0]     = true;
	component_mask[1]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  2,
						  ZeroFunction<dim>(dim+1+dim),
						  boundary_values,
						  component_mask);

	component_mask[0]     = false;
	component_mask[1]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  3,
						  ZeroFunction<dim>(dim+1+dim),
						  boundary_values,
						  component_mask);

	component_mask[0]     = false;
	component_mask[1]     = false;
	VectorTools::interpolate_boundary_values (dof_handler,
						  4,
						  ZeroFunction<dim>(dim+1+dim),
						  boundary_values,
						  component_mask);

	component_mask[0]     = true;
	component_mask[1]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  5,
						  NonhomDirichletBoundaryValues<dim>(time, test_case,alpha_eps),
						  boundary_values,
						  component_mask);

	component_mask[0]     = false;
	component_mask[1]     = false;
	component_mask[2]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  1,
						  //ZeroFunction<dim>(dim+1),
						  ConstantFunction<dim>(1.0,dim+1+dim),
						  boundary_values,
						  component_mask);


	VectorTools::interpolate_boundary_values (dof_handler,
						  5,
						  ConstantFunction<dim>(1.0,dim+1+dim),
						  boundary_values,
						  component_mask);

      }
    else if (test_case == "pressurized" || test_case == "Sneddon")
      {
	component_mask[0]     = true;
	component_mask[1]     = true;
	VectorTools::interpolate_boundary_values (dof_handler,
						  0,
						  ZeroFunction<dim>(dim+1+dim),
						  boundary_values,
						  component_mask);
      }
    else if (test_case == "screw_domi")
      {
	    component_mask[0]     = false;
    component_mask[1]     = false;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      ZeroFunction<dim>(dim+1+dim),
					      boundary_values,
					      component_mask);  

 
    component_mask[0]     = false;
    component_mask[1]     = false;
    VectorTools::interpolate_boundary_values (dof_handler,
					      1,
					      ZeroFunction<dim>(dim+1+dim),					
					      boundary_values,
					      component_mask);
 

    component_mask[0]     = false; 
    component_mask[1]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
					      ZeroFunction<dim>(dim+1+dim), 	
                                              boundary_values,
                                              component_mask);
 
    // TODO
    component_mask[0]     = false;  // false
    component_mask[1]     = true; 
    VectorTools::interpolate_boundary_values (dof_handler,
					      3,
					      NonhomDirichletBoundaryValues<dim>(time, test_case,alpha_eps),
					      boundary_values,
					      component_mask);

    component_mask[0]     = false;
    component_mask[1]     = false;
    component_mask[2]     = false;
    VectorTools::interpolate_boundary_values (dof_handler,
					      2,
					      //ZeroFunction<dim>(dim+1),
					      ConstantFunction<dim>(1.0,dim+1+dim),
					      boundary_values,
					      component_mask);


      }
    else if ((test_case == "Sneddon3D") || (test_case == "Het3D"))
      {
    component_mask[0]     = true;
    component_mask[1]     = true;
    component_mask[2]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      ZeroFunction<dim>(dim+1+dim),
					      boundary_values,
					      component_mask);  

 
    component_mask[0]     = true;
    component_mask[1]     = true;
    component_mask[2]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
					      1,
					      ZeroFunction<dim>(dim+1+dim),					
					      boundary_values,
					      component_mask);
 

    component_mask[0]     = true; 
    component_mask[1]     = true;
    component_mask[2]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
					      ZeroFunction<dim>(dim+1+dim), 	
                                              boundary_values,
                                              component_mask);
 
    component_mask[0]     = true; 
    component_mask[1]     = true; 
    component_mask[2]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
					      3,
					      ZeroFunction<dim>(dim+1+dim),  
					      boundary_values,
					      component_mask);

   component_mask[0]     = true; 
    component_mask[1]     = true; 
    component_mask[2]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
					      4,
					      ZeroFunction<dim>(dim+1+dim),  
					      boundary_values,
					      component_mask);

   component_mask[0]     = true; 
    component_mask[1]     = true; 
    component_mask[2]     = true;
    VectorTools::interpolate_boundary_values (dof_handler,
					      5,
					      ZeroFunction<dim>(dim+1+dim),  
					      boundary_values,
					      component_mask);
	

      }

    
    for (typename std::map<unsigned int, double>::const_iterator
	   i = boundary_values.begin();
	 i != boundary_values.end();
	 ++i)
      solution(i->first) = i->second;
    
}

// This function applies boundary conditions 
// to the Newton iteration steps. For all variables that
// have Dirichlet conditions on some (or all) parts
// of the outer boundary, we apply zero-Dirichlet
// conditions, now. 
template <int dim>
void
Dynamic_Fracture_Problem<dim>::set_newton_bc ()
{
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
	component_mask[1]     = false;
	component_mask[dim+1]     = true;
	// changed component_mask[dim+2]     = true;
  component_mask[dim+2]     = false;
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

// In this function, we solve the linear systems
// inside the nonlinear Newton iteration. We only
// use a direct solver from UMFPACK.
template <int dim>
void 
Dynamic_Fracture_Problem<dim>::solve () 
{
  timer.enter_section("Solve linear system.");
  Vector<double> sol, rhs;    
  sol = newton_update;    
  rhs = system_rhs;
  
  SparseDirectUMFPACK A_direct;
  A_direct.factorize(system_matrix);     
  A_direct.vmult(sol,rhs); 
  newton_update = sol;
  
  constraints.distribute (newton_update);
  timer.exit_section();
}

// This is the Newton iteration to solve the 
// non-linear system of equations. First, we declare some
// standard parameters of the solution method. Addionally,
// we also implement an easy line search algorithm. 
template <int dim>
double Dynamic_Fracture_Problem<dim>::newton_iteration (const double time) 
					       
{ 
  Timer timer_newton;
  //const double lower_bound_newton_residuum = 1.0e-8; 
  //const unsigned int max_no_newton_steps  = 20;

  // Decision whether the system matrix should be build
  // at each Newton step
  const double nonlinear_rho = 0.1; 
 
  // Line search parameters
  unsigned int line_search_step;

  // For residual-based Newton with increasing
  // residual use a small number of steps: 4
  //const unsigned int  max_no_line_search_steps = 10;
  const double line_search_damping = 0.6;
  double new_newton_residuum;
  
  // Application of the initial boundary conditions to the 
  // variational equations:
  // Removed because of coupled Newton iteration.
  assemble_system_rhs();

  double newton_residuum = system_rhs.linfty_norm(); 
  double old_newton_residuum= newton_residuum;
  unsigned int newton_step = 1;

  if (bool_set_explicitely_delta_fp)
    delta_fixed_point_newton = 0.0;
   
  if (newton_residuum < lower_bound_newton_residuum)
    {
      std::cout << '\t' 
		<< std::scientific 
		<< newton_residuum 
		<< std::endl;     
    }
  
  while (newton_residuum > lower_bound_newton_residuum &&
	 newton_step < max_no_newton_steps)
    {
      timer_newton.start();

  if (bool_set_explicitely_constant_k)
    {
      //delta_fixed_point_newton = 0.0;
      //delta_fixed_point_newton = newton_residuum;
      constant_k = newton_residuum;
      if (constant_k > (0.5 * alpha_eps))
	constant_k = 0.5 * alpha_eps;
      //std::cout << constant_k << std::endl;
    }
 

      old_newton_residuum = newton_residuum;
      
      assemble_system_rhs();
      newton_residuum = system_rhs.linfty_norm();

      if (newton_residuum < lower_bound_newton_residuum)
	{
	  std::cout << '\t' 
		    << std::scientific 
		    << newton_residuum << std::endl;
	  break;
	}

      if (newton_residuum > 1e+14)
	{
	  std::cout << '\t' 
		    << std::scientific 
		    << newton_residuum << std::endl;
	  std::cout << "Newton residual too high. Aborting." << std::endl;
	  abort();
	}


      // TODO: This seems to lead to 
      // non-robust Newton iterations
      //if (newton_residuum/old_newton_residuum > nonlinear_rho)
	assemble_system_matrix ();	

      // Solve Ax = b
      solve ();	  
        
      line_search_step = 0;	  
      for ( ; 
	    line_search_step < max_no_line_search_steps; 
	    ++line_search_step)
	{	 
	  solution +=newton_update;

	  assemble_system_rhs ();			
	  new_newton_residuum = system_rhs.linfty_norm();
	  
	  if (new_newton_residuum < newton_residuum)
	      break;
	  else 	  
	    solution -= newton_update;
	  
	  newton_update *= line_search_damping;

	}

      // Allow for increasing residual
      if (line_search_step == max_no_line_search_steps)
	{
	  solution +=newton_update;	   

	}
     
      timer_newton.stop();
      
      std::cout << std::setprecision(5) <<newton_step << '\t' 
		<< std::scientific << newton_residuum << '\t'
		<< std::scientific << newton_residuum/old_newton_residuum  <<'\t' ;
      if (newton_residuum/old_newton_residuum > nonlinear_rho)
	std::cout << "r" << '\t' ;
      else 
	std::cout << " " << '\t' ;
      std::cout << line_search_step  << '\t' 
		<< delta_fixed_point_newton  << '\t'
		<< std::scientific << timer_newton.wall_time ()
		<< std::endl;

//      if ((newton_residuum/old_newton_residuum > upper_newton_rho) && (newton_step > 1))
//	{
//	  break;
//	}


      if (bool_use_modified_Newton)
	{
	  // Update delta for dynamic switch between fixed point and Newton
	  double Qn = newton_residuum/old_newton_residuum;
	  double Qn_inv = old_newton_residuum/newton_residuum;
	  
	  // Mandel's et al. formula (31)
	  //delta_fixed_point_newton = delta_fixed_point_newton * (0.2 + 4.0/(0.7 + std::exp(1.5 * Qn)));
	  
	  delta_fixed_point_newton = delta_fixed_point_newton * (a_fp/(std::exp(Qn_inv)) + b_fp/(std::exp(Qn)));

	  // Normalize delta
	  if (delta_fixed_point_newton > 1.0)
	    delta_fixed_point_newton = 1.0;
	  else if (delta_fixed_point_newton < 0.0)
	    delta_fixed_point_newton = 0.0;
	  
	}

     

      // Updates
      timer_newton.reset();
      newton_step++;      
    }

  //std::cout << "NumNewtonIter: " << newton_step << std::endl;
  total_number_newton_steps += newton_step;

  return newton_residuum/old_newton_residuum; 
}


template <int dim>
void Dynamic_Fracture_Problem<dim>::newton_iteration_error_based (const double time) 
					       
{ 
   Timer timer_newton;
   //const double lower_bound_newton_update = 1.0e-8; //1.0e-12; 
   //const double lower_bound_newton_residuum = 1.0e-8; 
   //const unsigned int max_no_newton_steps  = 50;

  const double lambda_newton_lower_bound = 1.0e-10;

  // Decision whether the system matrix should be build
  // at each Newton step
  const double nonlinear_rho = 0.1; 
 
  // Line search parameters
  unsigned int line_search_step;
  const unsigned int  max_no_line_search_steps = 10;
  const double line_search_damping = 0.6;
  double new_newton_residuum;

  if (bool_set_explicitely_delta_fp)
    delta_fixed_point_newton = 0.0;

  
  // Application of the initial boundary conditions to the 
  // variational equations:
  //set_initial_bc (time);
  assemble_system_rhs();

  double newton_residuum = system_rhs.l2_norm(); 
  double newton_update_norm = 1.0;
  double old_newton_update_norm = newton_update_norm;
  double bar_newton_update_norm = newton_update_norm;
  double old_bar_newton_update_norm = newton_update_norm;
  unsigned int newton_step = 1;

  BlockVector<double> bar_newton_update = newton_update;
  BlockVector<double> old_bar_newton_update = bar_newton_update;
  BlockVector<double> new_newton_update = newton_update;
  double theta_newton = 1.0;
  double old_theta_newton = theta_newton;

  double lambda_newton = 1.0;
  double lambda_newton_prime = lambda_newton;
  double old_lambda_newton = lambda_newton;
  double tmp_lambda_newton = lambda_newton;

  double mu_newton = 1.0;
  double mu_newton_prime = 1.0;

  while (newton_update_norm > lower_bound_newton_update &&
	 newton_residuum > lower_bound_newton_residuum &&
	 newton_step < max_no_newton_steps)
    {
      timer_newton.start();

      // Step 1
      old_newton_update_norm = newton_update_norm;
      old_theta_newton = theta_newton;

      old_bar_newton_update = bar_newton_update;
      old_bar_newton_update_norm = bar_newton_update_norm;

      old_lambda_newton = lambda_newton;

      
      assemble_system_rhs();
 
      newton_residuum = system_rhs.l2_norm();

      if (newton_residuum < lower_bound_newton_residuum)
	{
	  std::cout << '\t' 
		    << std::scientific 
		    << newton_residuum << std::endl;
	  break;
	}





      //if (theta_newton > nonlinear_rho)
	assemble_system_matrix ();	


      // Solve Ax = b
      solve ();	

      new_newton_update = newton_update;
      newton_update_norm = new_newton_update.l2_norm(); 
      //      std::cout << newton_update_norm << std::endl;

      if (newton_update_norm < lower_bound_newton_update)
	{
	  solution += new_newton_update;
	  std::cout << '\t' 
		    << std::scientific 
		    << newton_update_norm << std::endl;
	  break;
	}

      

      BlockVector<double> diff_newton = bar_newton_update;
      for (unsigned int i=0; i < new_newton_update.size();i++)
	diff_newton[i] = old_bar_newton_update[i] - new_newton_update[i];

      double diff_newton_norm = diff_newton.l2_norm();

      mu_newton = 1.0;
      if (newton_step > 1)
      	mu_newton = (old_newton_update_norm * old_bar_newton_update_norm) / 
      	  (diff_newton_norm * newton_update_norm) * old_lambda_newton;
      

      lambda_newton = std::min(1.0, mu_newton);


    regu_test:
      //std::cout << "lambda newton: " << lambda_newton << std::endl;
      if (lambda_newton < lambda_newton_lower_bound)
	{
	  std::cout << lambda_newton << "\tConvergence failure. Aborting in Step 1." << std::endl;

	  // Option 1
	  //abort();	  

	  /*
	  // Option 2: switch to residual-based Newton
	  std::cout << lambda_newton << "\tConvergence failure in error-based Newton. Go to residual-based Newton." << std::endl;

	  // TODO: be careful, in standard newton the 
	  // initial bc are commented!!!!!
	  newton_iteration (time);
	  break;
	  */	  

	  
	  // Option 3: just take the small step and continue
	  std::cout << lambda_newton << "\tConvergence failure in error-based Newton. Take the smallest update and continue." << std::endl;
	  for (unsigned int i=0; i < new_newton_update.size();i++)
	    {
	      solution(i) = solution(i) + lambda_newton * new_newton_update(i);
	    }
	  break;	
	  

	}

      // Step 2
    step_2:

      tmp_lambda_newton = lambda_newton;
      for (unsigned int i=0; i < new_newton_update.size();i++)
	{
	  solution(i) = solution(i) + lambda_newton * new_newton_update(i);
	}

      assemble_system_rhs ();

      // Solve now the bar-system
      // old Jacobian and new right hand side
      solve ();
      bar_newton_update = newton_update;
      bar_newton_update_norm = bar_newton_update.l2_norm(); 
   
      // Step 3
      // Compute monitoring quantities
      BlockVector<double> diff_mu_newton_prime = new_newton_update;
      for (unsigned int i=0; i < new_newton_update.size();i++)
	{
	  diff_mu_newton_prime(i) =
	    bar_newton_update(i) - (1.0 - lambda_newton) * new_newton_update(i);
	}


      double diff_mu_newton_prime_norm = diff_mu_newton_prime.l2_norm();
 
      theta_newton = bar_newton_update_norm / newton_update_norm;
      //std::cout << "theta_newton: " <<  theta_newton << std::endl;


      mu_newton_prime = 1.0;
      if (newton_step > 1)
	mu_newton_prime = 0.5 * (newton_update_norm * lambda_newton * lambda_newton) /
	  diff_mu_newton_prime_norm;

      lambda_newton_prime = lambda_newton;
      double eps_theta_newton = 1.0; // 1.001
      if (theta_newton >= eps_theta_newton) //TODO 1.0
       {
	 lambda_newton_prime = std::min(mu_newton_prime, 0.5 * lambda_newton);
	 lambda_newton = lambda_newton_prime;
	 //std::cout << "Monitoring test: " << lambda_newton << std::endl;

	 for (unsigned int i=0; i < new_newton_update.size();i++)
	   {
	     solution(i) = solution(i) - tmp_lambda_newton * new_newton_update(i);
	   }
	 

	 goto regu_test;
	 
       }
     else 
       lambda_newton_prime = std::min(mu_newton_prime, 1.0);
        

      
     if ((lambda_newton < (lambda_newton_prime + 1.0e-10)) && 
	 (lambda_newton > (lambda_newton_prime - 1.0e-10)) && 
	 (lambda_newton < 1.0 + 1.0e-10) &&
	 (lambda_newton > 1.0 - 1.0e-10)
	 )
       {
//	 std::cout << "Choice 1: " 
//		   << lambda_newton << "  " 
//		   << bar_newton_update_norm << "  "
//		   << theta_newton << std::endl;

	 if (bar_newton_update_norm < lower_bound_newton_update)
	   {
	     solution += bar_newton_update;
	     std::cout << '\t' 
		       << std::scientific 
		       << bar_newton_update_norm << std::endl;
	     break;
	   }
	 
       }
//     else if (lambda_newton_prime >= 4.0 * lambda_newton)
//       {
//	 // std::cout << "Choice 2:" << lambda_newton_prime << "   "  << lambda_newton  << std::endl;
//
//	 for (unsigned int i=0; i < new_newton_update.size();i++)
//	   {
//	     solution(i) = solution(i) - tmp_lambda_newton * new_newton_update(i);
//	   }
//	 
//	 
//	 lambda_newton = lambda_newton_prime;
//
//	 goto step_2;
//       }
        
     
     //std::cout << "Choice 3: " << lambda_newton << std::endl;

      assemble_system_rhs();
      newton_residuum = system_rhs.l2_norm();

      timer_newton.stop();
      
      std::cout << std::setprecision(5) <<newton_step << '\t' 
		<< std::scientific << newton_update_norm << '\t'
		<< std::scientific << bar_newton_update_norm << '\t'
		<< std::scientific << newton_residuum << '\t'
		<< std::scientific << newton_update_norm/old_newton_update_norm  <<'\t' ;
      if (old_theta_newton > nonlinear_rho)
	std::cout << "r" << '\t' ;
      else 
	std::cout << " " << '\t' ;
        
      std::cout << delta_fixed_point_newton  << '\t'
		<< std::scientific << timer_newton.wall_time ()
		<< std::endl;


	if (bool_use_modified_Newton)
	{
	  // Update delta for dynamic switch between fixed point and Newton
	  double Qn = bar_newton_update_norm/newton_update_norm;
	  double Qn_inv = newton_update_norm/bar_newton_update_norm;
	  
	  // Mandel's et al. formula (31)
	  //delta_fixed_point_newton = delta_fixed_point_newton * (0.2 + 4.0/(0.7 + std::exp(1.5 * Qn)));
	  
	  delta_fixed_point_newton = delta_fixed_point_newton * (a_fp/(std::exp(Qn_inv)) + b_fp/(std::exp(Qn)));

	  // Normalize delta
	  if (delta_fixed_point_newton > 1.0)
	    delta_fixed_point_newton = 1.0;
	  else if (delta_fixed_point_newton < 0.0)
	    delta_fixed_point_newton = 0.0;
	  
	}




      // Updates
      timer_newton.reset();
      newton_step++;      
    }

  //std::cout << "NumNewtonIter: " << newton_step << std::endl;
  total_number_newton_steps += newton_step;

}




template <int dim>
void
Dynamic_Fracture_Problem<dim>::project_back_phase_field ()  
{
  for (unsigned int i=0; i<solution.block(1).size(); ++i)
    if (solution.block(1)(i) < 0)
      solution.block(1)(i) = 0;
    else if (solution.block(1)(i) > 1)
      solution.block(1)(i) = 1;
}



// This function is known from almost all other 
// tutorial steps in deal.II.
template <int dim>
void
Dynamic_Fracture_Problem<dim>::output_results (const unsigned int refinement_cycle,
			      const BlockVector<double> output_vector)  const
{

  std::vector<std::string> solution_names;
  solution_names.push_back ("x_dis");
  solution_names.push_back ("y_dis");
  if (dim == 3)
    solution_names.push_back ("z_dis");
  solution_names.push_back ("phase_field");

  solution_names.push_back ("x_vel");
  solution_names.push_back ("y_vel");
  if (dim == 3)
    solution_names.push_back ("z_vel");

 std::vector<std::string> solution_names_lambda_penal;
  solution_names_lambda_penal.push_back ("x_lp");
  solution_names_lambda_penal.push_back ("y_lp");
  if (dim == 3)
    solution_names_lambda_penal.push_back ("z_lp");

  solution_names_lambda_penal.push_back ("pf_lp");

  solution_names_lambda_penal.push_back ("vx_lp");
  solution_names_lambda_penal.push_back ("vy_lp");
  if (dim == 3)
    solution_names_lambda_penal.push_back ("vz_lp");
      

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim+1+dim, DataComponentInterpretation::component_is_scalar);


  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);  
   
  data_out.add_data_vector (output_vector, solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  
  // only required for heterogeneous materials
  //data_out.add_data_vector(lame_coefficient_mu_vector, "mu");
  //data_out.add_data_vector(lame_coefficient_lambda_vector, "lambda");

  if (test_case == "screw_domi")
    data_out.add_data_vector(solution_stress_per_cell, "stress"); 


  data_out.add_data_vector (solution_lambda_penal_func, solution_names_lambda_penal,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);


  Vector<float> e_mod(triangulation.n_active_cells());
 if (test_case == "Het3D")
    {
 
      typename DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();

      unsigned int cellindex = 0;
      for (; cell != endc; ++cell, ++cellindex)
          {
            e_mod(cellindex) = 1.0 + func_emodulus->value(cell->center(), 0);
          }
      data_out.add_data_vector(e_mod, "E_modulus");
    }


  data_out.build_patches ();

  // Output
  // 1X XXXX := predictor output
  // 2X XXXX := corrector output
  // 3X XXXX := old solution on new mesh
  
  // X0 XXXX := level 0 of pred-corr
  // X1 XXXX := level 1 of pred-corr
  // etc.
   
  std::ostringstream filename;

  std::cout << "------------------" << std::endl;
  std::cout << "Write solution" << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << std::endl;
  filename << filename_basis
	   << Utilities::int_to_string (refinement_cycle, 6)
	   << ".vtk";
  
  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);

}

// With help of this function, we extract 
// point values for a certain component from our
// discrete solution. We use it to gain the 
// displacements of the structure in the x- and y-directions.
template <int dim>
double Dynamic_Fracture_Problem<dim>::compute_point_value (Point<dim> p, 
					       const unsigned int component) const  
{
 
  Vector<double> tmp_vector(dim+1);
  VectorTools::point_value (dof_handler, 
			    solution, 
			    p, 
			    tmp_vector);
  
  return tmp_vector(component);
}

// Now, we arrive at the function that is responsible 
// to compute the line integrals for stress evaluations
template <int dim>
void Dynamic_Fracture_Problem<dim>::compute_functional_values()
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
 if (test_case == "miehe_shear")
   bc_color = 3;
 else if (test_case == "l_shaped")
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


   if (test_case == "miehe_shear" || 
       test_case == "miehe_tension")
     {
       load_value[0] *= -1.0;

       std::cout << "Stress(x): " << time << "   " << load_value[0] << std::endl;
       std::cout << "Stress(y): " << time << "   " << load_value[1] << std::endl;

     }
   else if (test_case == "l_shaped")
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
   
       std::cout << "Stress(x): " << time << "\t" << load_increment << "\t" << load_value[0] << std::endl;
       std::cout << "Stress(y): " << time << "\t" << load_increment << "\t" << load_value[1] << std::endl;

     }
}



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
 if (test_case == "miehe_shear")
   bc_color = 3;
 else if (test_case == "l_shaped")
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


   if (test_case == "miehe_shear" || 
       test_case == "miehe_tension")
     {
       load_value[0] *= -1.0;


     }
   else if (test_case == "l_shaped")
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

template <int dim>
void Dynamic_Fracture_Problem<dim>::compute_stress_per_cell ()
{

  QGauss<dim> quad_gauss(degree+2);
  FEValues<dim> fe_values (fe, quad_gauss, update_gradients | update_quadrature_points);			   
  const unsigned int  n_q_points  = quad_gauss.size();
 std::vector<std::vector<Tensor<1,dim> > > old_solution_grads (n_q_points, 
								std::vector<Tensor<1,dim> > (dim+1+dim));
  




 typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

 //unsigned int cell_counter;

  solution_stress_per_cell.reinit(triangulation.n_active_cells());

  // This is the identity matrix in two dimensions:
  const Tensor<2,dim> Identity = ALE_Transformations
    ::get_Identity<dim> ();
 

  for (unsigned int cell_counter = 0; cell!=endc; ++cell, ++cell_counter)
    {
      fe_values.reinit (cell);
      fe_values.get_function_gradients (solution, old_solution_grads);

      double norm_stress_term = 0.0;
      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  const Tensor<2,dim> grad_u = ALE_Transformations
	    ::get_grad_u<dim> (q, old_solution_grads);

	   const Tensor<2,dim> E = 0.5 * (grad_u + transpose(grad_u));
	   const double tr_E = grad_u[0][0] + grad_u[1][1];


	  Tensor<2,dim> stress_term;
	  stress_term.clear();
	  stress_term = lame_coefficient_lambda * tr_E * Identity +  2 * lame_coefficient_mu * E;

	  norm_stress_term += ALE_Transformations
	    ::get_deviator_norm<dim> (stress_term);



	}

      solution_stress_per_cell(cell_counter) = (norm_stress_term  / n_q_points);
      

    }

}

template <int dim>
void
Dynamic_Fracture_Problem<dim>::compute_energy()
{
  // What are we computing? In Latex-style it is:
  // bulk energy = [(1+k)phi^2 + k] psi(e)
  // crack energy = \frac{G_c}{2}\int_{\Omega}\Bigl( \frac{(\varphi - 1)^2}{\eps}
  //+ \eps |\nabla \varphi|^2 \Bigr) \, dx

  const QGauss<dim> quadrature_formula(degree+2);
  const unsigned int n_q_points = quadrature_formula.size();

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_quadrature_points | update_JxW_values
                          | update_gradients);

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active(), endc = dof_handler.end();

 std::vector<Vector<double> > 
    old_solution_values (n_q_points, Vector<double>(dim+1+dim));

  std::vector<std::vector<Tensor<1,dim> > > 
    old_solution_grads (n_q_points, std::vector<Tensor<1,dim> > (dim+1+dim));


  double local_bulk_energy = 0.0;
  double local_crack_energy = 0.0;

  for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);


        fe_values.get_function_values(solution, old_solution_values);
        fe_values.get_function_gradients(solution, old_solution_grads);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
	    const Tensor<2,dim> grad_u = ALE_Transformations 
	      ::get_grad_u<dim> (q, old_solution_grads);

            const Tensor<2,dim> E = 0.5 * (grad_u + transpose(grad_u));
            const double tr_E = dealii::trace(E);

	    const double pf = old_solution_values[q](dim);

            const double tr_e_2 = dealii::trace(E*E);

            const double psi_e = 0.5 * lame_coefficient_lambda * tr_E*tr_E + lame_coefficient_mu * tr_e_2;

            local_bulk_energy += ((1+constant_k)*pf*pf+constant_k) * psi_e * fe_values.JxW(q);

            local_crack_energy += G_c/2.0 * ((pf-1) * (pf-1)/alpha_eps + alpha_eps * dealii::scalar_product(grad_u, grad_u))
                                  * fe_values.JxW(q);
          }

      }

  double bulk_energy = local_bulk_energy;
  double crack_energy = local_crack_energy;;

  std::cout << "Energies: " << timestep_number << "  " << time
        << "   " << bulk_energy
        << "   " << crack_energy
	    << std::endl;

}



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
  else if (test_case == "Sneddon3D")
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
  else if (test_case == "Het3D")
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


  template <int dim>
  void Dynamic_Fracture_Problem<dim>::update_quadrature_point_history ()
  {
    QGauss<dim>   quadrature_formula(degree+2); 

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values | update_gradients | update_quadrature_points);
 
   std::vector<std::vector<Tensor<1,dim> > >
     old_solution_grads (quadrature_formula.size(),
			 std::vector<Tensor<1,dim> >(dim+1+dim));


   const Tensor<2, dim> Identity =
     ALE_Transformations::get_Identity<dim>();

   Tensor<2,dim> zero_matrix;
   zero_matrix.clear();



    for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
         cell != dof_handler.end(); ++cell)
      {
          PointHistory<dim> *local_quadrature_points_history
            = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
          Assert (local_quadrature_points_history >=
                  &quadrature_point_history.front(),
                  ExcInternalError());
          Assert (local_quadrature_points_history <
                  &quadrature_point_history.back(),
                  ExcInternalError());

          fe_values.reinit (cell);
          fe_values.get_function_gradients (solution,
                                            old_solution_grads);

          // Then loop over the quadrature points of this cell:
          for (unsigned int q=0; q<quadrature_formula.size(); ++q)
            {
              // On each quadrature point, compute the strain increment from
              // the gradients, and multiply it by the stress-strain tensor to
              // get the stress update. Then add this update to the already
              // existing strain at this point.

	     

	      const Tensor<2,dim> grad_u = ALE_Transformations 
		::get_grad_u<dim> (q, old_solution_grads);
	      

	      //std::cout << "Bin drin" << std::endl;
	      // Linearized strain
	      const Tensor<2,dim> E = 0.5 * (grad_u + transpose(grad_u));
	      
	      const double tr_E = Structure_Terms_in_ALE
		::get_tr_E<dim> (E);
	      
	      
	      Tensor<2,dim> stress_term;
	      stress_term.clear();
	      stress_term = lame_coefficient_lambda * tr_E * Identity
		+ 2 * lame_coefficient_mu * E;

	      Tensor<2,dim> stress_term_plus;
              Tensor<2,dim> stress_term_minus;
	      
	      if (timestep_number > 1)
		{
		  decompose_stress(stress_term_plus, stress_term_minus,
                                   E, tr_E, zero_matrix , 0.0,
                                   lame_coefficient_lambda,
                                   lame_coefficient_mu, false);
		}
	      else 
		{
		  stress_term_plus = stress_term;
		  stress_term_minus = 0;
		}
	      

	      if (bool_set_initial_strain_history && (test_case == "pressurized"))
		{
		  //std::cout << "Drin" << std::endl;
		  double top = 2.0 + alpha_eps/4.0;
		  double bottom = 2.0 - alpha_eps/4.0;
		  if (((fe_values.quadrature_point(q)[0] >= 1.8) && (fe_values.quadrature_point(q)[0] <= 2.2)) &&
		      ((fe_values.quadrature_point(q)[1] >= bottom) && (fe_values.quadrature_point(q)[1] <= top))
		      )
		    {
		      stress_term_plus[0][0] = 1.0e+10;
		      stress_term_plus[0][1] = 0.0;
		      stress_term_plus[1][0] = 0.0;
		      stress_term_plus[1][1] = 1.0e+10;
		    }
		  else 
		    {
		      stress_term_plus[0][0] = 0.0;
		      stress_term_plus[0][1] = 0.0;
		      stress_term_plus[1][0] = 0.0;
		      stress_term_plus[1][1] = 0.0;
		    }
		  
		}


	      if (bool_set_initial_strain_history && (test_case == "miehe_shear"))
		{
		  //std::cout << "Drin" << std::endl;
		   double top = 0.5 + alpha_eps/4.0;
		   double bottom = 0.5 - alpha_eps/4.0;
		 
		  if (((fe_values.quadrature_point(q)[0] >= 0.5) && (fe_values.quadrature_point(q)[0] <= 1.0)) &&
		      ((fe_values.quadrature_point(q)[1] >= bottom) && (fe_values.quadrature_point(q)[1]  <= top))
		      )
		    {
		      stress_term_plus[0][0] = 0.0; //1.0e+14;
		      stress_term_plus[0][1] = 0.0;
		      stress_term_plus[1][0] = 0.0;
		      stress_term_plus[1][1] = 0.0; //1.0e+14;
		    }
		  else 
		    {
		      stress_term_plus[0][0] = 0.0;
		      stress_term_plus[0][1] = 0.0;
		      stress_term_plus[1][0] = 0.0;
		      stress_term_plus[1][1] = 0.0;
		    }
		  
		}
	      	      


 
	      double stress_term_plus_norm = std::sqrt(stress_term_plus[0][0] * stress_term_plus[0][0] +
						       stress_term_plus[0][1] * stress_term_plus[0][1] +
						       stress_term_plus[1][0] * stress_term_plus[1][0] +
						       stress_term_plus[1][1] * stress_term_plus[1][1]);

	      Tensor<2,dim> old_stress = local_quadrature_points_history[q].old_stress;

	      double old_stress_norm = std::sqrt(old_stress[0][0] * old_stress[0][0] +
						 old_stress[0][1] * old_stress[0][1] +
						 old_stress[1][0] * old_stress[1][0] +
						 old_stress[1][1] * old_stress[1][1]);

	      Tensor<2,dim> new_stress;
	      if (stress_term_plus_norm > old_stress_norm)
		{
		  new_stress  = stress_term_plus; //old_stress + stress_term_plus;
		  //std::cout << "New stress" << std::endl;
		}
	      else
		{
		  new_stress  = old_stress;
		  //std::cout << "Oooooooooooooooold stress" << std::endl;
		}

              // The result of all these operations is then written back into
              // the original place:
              local_quadrature_points_history[q].old_stress
                = new_stress;
            }
        }
  }




template <int dim>
void Dynamic_Fracture_Problem<dim>::solve_spatial_problem () 
{ 

  BlockVector<double> tmp_old_sol;
  tmp_old_sol = solution;

  unsigned int prec_corr_counter_output_a = 100000 + timestep_number;
  unsigned int prec_corr_counter_output_b = 200000 + timestep_number;
  unsigned int prec_corr_counter_output_c = 300000 + timestep_number;

    redo_step:
      std::cout << "Timestep " << timestep_number 
		<< " (" << time_stepping_scheme 
		<< ")" <<    ": " << time
		<< " (" << timestep << ")"
		<< "   " << "Cells: " << triangulation.n_active_cells()
		<< "   " << "DoFs: " << dof_handler.n_dofs()
		<< "\n======================================" 
		<< "==========================================" 
		<< std::endl;
      if (!bool_use_error_oriented_Newton)
	{
	  std::cout << "Iter\t" << "Residual\t"
		    << "Red. New/Old\t" 
		    << "Bui.Mat\t"
		    << "LS\t"
		    << "Delta (Fix-P)\t"
		    << "CPU time"
		    << std::endl;
	}
      else if (bool_use_error_oriented_Newton)
	{
	  std::cout << "Iter\t" << "Update\t\t"
		    << "Update(bar)\t" 
		    << "Residual\t"
		    << "Red. Up/Old-Up\t"
		    << "Bui.Mat\t"
		    << "CPU time"
		    << std::endl;


	}


      BlockVector<double> solution_lambda_penal_func_difference;
      penal_iterations = 0;
      double L2_error = 0.0;
      double L2_error_relative = 0.0;
      double initial_L2_norm = 0.0;

      /*
	// TODO
      if (test_case == "miehe_shear" ||
	  test_case == "pressurized" ||
	  test_case == "screw_domi" ||
	  test_case == "Sneddon")
	{
	  lower_bound_newton_update = 1.0e-8;
	  lower_bound_newton_residuum = 1.0e-8; 
	}
      else if (test_case == "l_shaped")
	{
	  lower_bound_newton_update = 1.0e-8;
	  lower_bound_newton_residuum = 1.0e-8; 
	}
      */

      double pressure_increment_factor = 0.1;
      
 pf_extra_step:
      
      // Reset penal function for new time step
      solution_lambda_penal_func = 0;

      // Augmented Lagrangian iteration
      do
	{
	  // Step 1: Solve forward system
	  set_initial_bc (time);
	  if (!bool_use_error_oriented_Newton)
	    {
	      delta_fixed_point_newton = 1.0;
	      newton_iteration (time);
	    }
	  else if (bool_use_error_oriented_Newton)
	    {
	      // Only for initialization - apart from this,
	      // delta_fixed_point Newton is not used
	      // in the error oriented Newton method
	      delta_fixed_point_newton = 1.0;
	      newton_iteration_error_based (time);   
	    }
 

	  
	  // Step 2: Update lambda
	  old_timestep_solution_lambda_penal_func = solution_lambda_penal_func;
	  for (unsigned int j=0; j<solution_lambda_penal_func.size(); ++j)
	    {
	      if ((solution_lambda_penal_func(j)  + gamma_penal * (solution(j)  - old_timestep_solution(j) )) > 0.0)
		solution_lambda_penal_func(j) 
		  = solution_lambda_penal_func(j)  + gamma_penal * (solution(j)  - old_timestep_solution(j) );
	      else 
		solution_lambda_penal_func(j) = 0.0;

	      //std::cout << solution_lambda_penal_func(j) << std::endl;
	    }

	  // Step 3: Check stopping criterium
	  solution_lambda_penal_func_difference = solution_lambda_penal_func;
	  solution_lambda_penal_func_difference -= old_timestep_solution_lambda_penal_func;
	  

	  ComponentSelectFunction<dim> value_select (dim, dim+1+dim);
	  Vector<float> difference_per_cell (triangulation.n_active_cells());
	  VectorTools::integrate_difference (dof_handler,
					     //tmp_vec,
					     solution_lambda_penal_func_difference,
					     ZeroFunction<dim>(dim+1+dim),					     
					     difference_per_cell,
					     QGauss<dim>(fe.degree+1),
					     VectorTools::L2_norm,
					     &value_select);
	  L2_error = difference_per_cell.l2_norm();

	  if (penal_iterations == 0)
	    initial_L2_norm = L2_error;
	  
	  L2_error_relative = L2_error / initial_L2_norm;



	  if (bool_use_adaptive_newton_bound && (L2_error_relative > 1.0e-4))
	    {
	      lower_bound_newton_update = 1.0e-4 * L2_error_relative;
	      lower_bound_newton_residuum = 1.0e-4 * L2_error_relative; 
	    }
	  else
	    {
	      lower_bound_newton_update = 1.0e-8; 
	      lower_bound_newton_residuum = 1.0e-8; 

	      if (test_case == "miehe_shear")
		{
		  lower_bound_newton_update = 1.0e-8;
		  lower_bound_newton_residuum = 1.0e-8; 
		}

	    }



	  std::cout << std::endl;
	  std::cout << "AL error: " 
		    << penal_iterations << "\t" 
		    << solution_lambda_penal_func.l2_norm() << "\t" 
		    << L2_error << "\t"
		    << L2_error_relative << "\t"
		    << lower_bound_newton_residuum
		    << std::endl;
	  std::cout << std::endl; 

	  
	  penal_iterations++; 

	}
      while ((L2_error > tolerance_absolute_augmented_L_iteration) &&
	     (L2_error_relative > tolerance_augmented_L_iteration) && 
	     (penal_iterations < max_no_of_augmented_L_penal_iterations));
      //      std::cout << std::endl;
      std::cout << "NumALIter: " << penal_iterations << std::endl;



  
      /*  
	  // TODO: evtl wieder einkommentieren
      // Option 1 (check at each time step first for the larger pressure increment)
      if ((test_case == "pressurized") && 
	  (penal_iterations == max_no_of_augmented_L_penal_iterations))
	{

	  solution = tmp_old_sol;
	  penal_iterations = 0;
	  current_pressure -= pressure_increment;
	  current_pressure = current_pressure + pressure_increment_factor * pressure_increment;
	  pressure_increment_factor *= 0.1; // TODO 0.1
	  if (pressure_increment_factor < 1e-6)
	    {
	      std::cout << "Aborting. Too small pressure_increment." << std::endl;
	      abort();
	    }
	  
	  std::cout << "Repeat step with smaller pressure increment: " << time << "\t" << current_pressure << "\t" << pressure_increment_factor << std::endl;

	  goto pf_extra_step;
	}
      */


/*
      // Option 2 (once smaller increment is used, keep this smaller one)
      if ((test_case == "pressurized") && 
	  (penal_iterations == max_no_of_augmented_L_penal_iterations))
	{

	  solution = tmp_old_sol;
	  penal_iterations = 0;
	  current_pressure -= pressure_increment;
	  pressure_increment =  pressure_increment_factor * pressure_increment
	  current_pressure = current_pressure + pressure_increment;

	  pressure_increment_factor *= 0.1;

	  if (pressure_increment_factor < 1e-6)
	    {
	      std::cout << "Aborting. Too small pressure_increment." << std::endl;
	      abort();
	    }
	  
	  std::cout << "Repeat step with smaller pressure increment: " << time << "\t" << pressure_increment << "\t" << pressure_increment_factor << std::endl;

	  goto pf_extra_step;
	}
*/            



      // Nov 21, 2016: Ist bei predictor-corrector refinement
      // und augmented Lagrangian scheinbar hier notwendig
      // Entscheidend war aber eine hoehere Anzahl von line search
      // steps (anstatt 4 nehme ich wieder 10)
      project_back_phase_field();
          
      // TODO
      // Predictor-correcter check
      if (test_case == "miehe_shear" ||
	  test_case == "l_shaped")
	{
	  // Plot predictor solution on old mesh
	  if (bool_plot_additional_solutions)
	    {
	      output_results (prec_corr_counter_output_a,solution);
	      prec_corr_counter_output_a += 10000;
	    }

	  bool changed = refine_mesh();
	  if (changed)
	    {
	      // redo the current time step
	      std::cout << "Mesh change: corrector step" << std::endl;

	      // Plot predictor solution on new mesh
	      if (bool_plot_additional_solutions)
		{
		  output_results (prec_corr_counter_output_b,solution);
		  prec_corr_counter_output_b += 10000;
		}

	      //time -= timestep;
	      solution = old_timestep_solution;
	      
	      // Plot old solution on new mesh
	      if (bool_plot_additional_solutions)
		{
		  output_results (prec_corr_counter_output_c,solution);
		  prec_corr_counter_output_c += 10000;
		}

	      goto redo_step;

	    }
	}


	// Total number of Newton steps per time step
      std::cout << "NumNewtonIterTotal: " << total_number_newton_steps << std::endl;
      std::cout << std::endl;

}



// As usual, we have to call the run method. It handles
// the output stream to the terminal.
// Second, we define some output skip that is necessary 
// (and really useful) to avoid to much printing 
// of solutions. For large time dependent problems it is 
// sufficient to print only each tenth solution. 
// Third, we perform the time stepping scheme of 
// the solution process.
template <int dim>
void Dynamic_Fracture_Problem<dim>::run () 
{ 
  // Switch dimension !!

  // Defining test cases
  test_case = "miehe_shear";
  //test_case = "l_shaped";
  //test_case = "Sneddon";
  // before it was pressurized test case
  // test_case = "pressurized";
  //test_case = "screw_domi";
  //test_case = "Sneddon3D";
  //test_case = "Het3D";

  // Switch dimension !!
  

  // Initialization of some variables
  // that may be overwritten in the set_runtime routines
  set_global_parameters ();

  // Setting specific parameters for each problem
  if (test_case == "miehe_tension" ||
      test_case == "miehe_shear")
    set_runtime_parameters_Miehe ();
  else if (test_case == "l_shaped")
    set_runtime_parameters_L_shaped ();
  else if (test_case == "Sneddon")
    set_runtime_parameters_Sneddon ();
  else if (test_case == "pressurized")
    set_runtime_parameters_pressurized ();
  else if (test_case == "screw_domi")
    set_runtime_parameters_screw_domi ();
  else if (test_case == "Sneddon3D")
    {
      set_runtime_parameters_Sneddon3D ();
      // TODO: set dimension to 3 !!
    }
  else if (test_case == "Het3D")
    {
      set_runtime_parameters_Het3D ();
      // TODO: set dimension to 3 !!
    }
  else 
    {
      std::cout << "Framework not implemented. Aborting. " << std::endl;
      abort();
    }

  setup_system();


  if (test_case == "Sneddon3D" 
      || test_case == "Het3D"
      )
    {
      for (unsigned int i=0;i<3;i++)
	refine_mesh();
    }


  double min_cell_diameter = 1.0e+10;
  /*
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  for (; cell!=endc; ++cell)
    { 
      cell_diameter = cell->diameter();
      if (min_cell_diameter > cell_diameter)
	min_cell_diameter = cell_diameter;	
      
    }
  */
  
  min_cell_diameter = dof_handler.begin(0)->diameter();
  min_cell_diameter *= std::pow(2.,-1.0*( global_refinement_steps  +  pred_corr_levels ));
  std::cout << "Min cell dia: " << min_cell_diameter << std::endl;


  // Set phase-field parameters
  constant_k = 1.0e-12;// * min_cell_diameter;
  // TODO - adjust according to mesh levels
  alpha_eps =  2.0 * min_cell_diameter; 

  //if (test_case == "miehe_shear")
  //  alpha_eps =  0.0884;
  

  std::cout << "\n==============================" 
	    << "====================================="  << std::endl;
   std::cout << "Parameters\n" 
	    << "==========\n"
	     << "min h:             "   << min_cell_diameter << "\n"
	     << "k:                 "   <<  constant_k << "\n"  
	     << "eps:               "   <<  alpha_eps << "\n"
	     << "G_c:               "   <<  G_c << "\n"
	     << "delta penal:       "   <<  delta_penal << "\n"
	     << "Poisson nu:        "   <<  poisson_ratio_nu << "\n"
	     << "Lame mu:           "   <<  lame_coefficient_mu << "\n"
	     << "Lame lambda:       "   <<  lame_coefficient_lambda << "\n"
	     << std::endl;
   
   
   for (unsigned int i=0;i< pred_corr_levels;i++)
     {
       if (test_case == "miehe_shear" ||
	   test_case == "miehe_tension" ||
	   test_case == "screw_domi" ||
	   test_case == "l_shaped")
	 {
	   ConstraintMatrix constraints;
	   constraints.close();
	   
	   std::vector<bool> component_mask (dim+1+dim, true);
	   
	   if (sub_test_case == "hohlwalzung")
	     {
	       VectorTools::project (dof_handler,
				     constraints,
				     QGauss<dim>(degree+2),
				     InitialValuesHohlwalzung<dim>(alpha_eps),
				     solution
				     );
	     }
	   else 
	     {
	       VectorTools::project (dof_handler,
				     constraints,
				     QGauss<dim>(degree+2),
				     InitialValuesMiehe<dim>(alpha_eps,bool_initial_crack_via_phase_field),
				     solution
				     );
	     }
	   
	   compute_stress_per_cell();
	   output_results (0,solution);
	 }
       else if (test_case == "pressurized" ||
		test_case == "Sneddon")
	 {
	   /*
	   ConstraintMatrix constraints;
	   constraints.close();
	   
	   std::vector<bool> component_mask (dim+1+dim, true);
	   VectorTools::project (dof_handler,
				 constraints,
				 QGauss<dim>(degree+2),
				 InitialValuesPressurized<dim>(alpha_eps),
				 solution
				 );
	   */
	   VectorTools::interpolate(dof_handler,
				    InitialValuesPressurized<dim>(alpha_eps), 
				    solution);

	   output_results (i,solution);
	 }
       else if (test_case == "Sneddon3D")
	 {
	   ConstraintMatrix constraints;
	   constraints.close();
	   
	   std::vector<bool> component_mask (dim+1+dim, true);
	   VectorTools::project (dof_handler,
				 constraints,
				 QGauss<dim>(degree+2),
				 InitialValuesSneddon3D<dim>(alpha_eps),
				 solution
				 );
	   
	   output_results (0,solution);
	 }
       else if (test_case == "Het3D")
	 {
	   ConstraintMatrix constraints;
	   constraints.close();
	   
	   std::vector<bool> component_mask (dim+1+dim, true);
	   VectorTools::project (dof_handler,
				 constraints,
				 QGauss<dim>(degree+2),
				 InitialValuesHet3D<dim>(alpha_eps),
				 solution
				 );
	   
	   output_results (0,solution);
	 }
       
       cout<<"PRE REFINEMENT STEP :: " << i << endl;
       refine_mesh();
      	 
       
     }

   // RESET THE INITIAL VALUES
   
   if (test_case == "pressurized" ||
       test_case == "Sneddon")
     {
        VectorTools::interpolate(dof_handler,
				    InitialValuesPressurized<dim>(alpha_eps), 
				    solution);

     }
   
   project_back_phase_field();

   
   
   if (bool_use_strain_history)
     {
       // Initialize strain history
       bool_set_initial_strain_history = true;
       update_quadrature_point_history ();
       output_results (0,solution);
       bool_set_initial_strain_history = false;
     }

  unsigned int refine_mesh_1st = 500000;
  unsigned int refine_mesh_2nd = 20000;
  unsigned int refine_mesh_3rd = 30000;
  unsigned int refine_mesh_4th = 40000;
  unsigned int refine_mesh_5th = 500000;
 
  const unsigned int output_skip = 1;
  unsigned int refinement_cycle = 0;

  double tmp_timestep;

  
  // Initialize old and old_old timestep sizes
  old_timestep = timestep;
  old_old_timestep = timestep;

  if (test_case == "pressurized")
    {
      current_pressure = 0.0; //1.0e-1; 
      old_timestep_current_pressure = 0.0; //1.0e-1; 
    }
  else 
    {
      current_pressure = 0.0;
      old_timestep_current_pressure = 0.0; //1.0e-1; 
    }

  double pressure_increment = 1.0e-1;

  double REL_functional_error = 0.0;
  double REL_old_timestep_functional_error = 0.0;



  // Time loop
  do
    { 

      // Possibly adapt time step size
      //if (timestep_number > 5)
      //	timestep = 1.0e-2;

      //     if (timestep_number > 6)
      //	timestep = 1.0e-8;
//      if (timestep_number > 9)
//	{
//	  pressure_increment = 0.01;
//	  //gamma_penal = 1.0e+4;
//	}
//      else
//	{
//	  //pressure_increment = 1.0e+3;
//	  //gamma_penal = 1.0e+5;
//	}


      old_timestep_current_pressure = current_pressure;
      if (test_case == "pressurized")
	current_pressure += pressure_increment; //5.0e-2 + time * 5.0e-2; //1.0e-1 + time * 1.0e-1;
      else if ((test_case == "Sneddon") || (test_case == "Sneddon3D"))
	current_pressure = 1.0e-3;
      else if (test_case == "Het3D")
	current_pressure = 1.0e-3 + time * 0.25;
      else 
	current_pressure = 0.0;

     

      if (test_case == "pressurized" || test_case == "Het3D")
	std::cout << "Current pressure: " << time << "\t" << current_pressure << std::endl;
      
      std::cout << std::endl;
      
      // Compute next time step
      total_number_newton_steps = 0;
      old_old_timestep_solution = old_timestep_solution;
      old_timestep_solution = solution;


      // Update time step sizes
      old_old_timestep = old_timestep;
      old_timestep = timestep;

      tmp_solution_for_time_adaptivity = solution;


      if (use_time_step_control)
	{
	  if (timestep_number <= number_of_nonadaptive_time_steps) 
	    time_stepping_scheme = "BE";
	  else
	    time_stepping_scheme = "BE_adaptive"; 

	  std::cout << "Not tested in detail." << std::endl;
	  abort();
	  
	}

 
     if (time_stepping_scheme == "BE")
	{
	  solve_spatial_problem ();
	  time += timestep;
	}
     else if (time_stepping_scheme == "CN")
       {
	 // n = 0 BE to obtain "smooth" solution
	 if (timestep_number == 0)
	   {
	     theta = 1.0;
	   }
	 else
	   theta = 0.5;
	
	 solve_spatial_problem ();
	 time += timestep;

	 
       }
      else if (time_stepping_scheme == "BE_adaptive")
      {
	tmp_timestep = timestep;
	double tmp_time = time;
	unsigned int counter_redo_timestep = 0;

	// Adaptive time step control loop
	do {
	  if (counter_redo_timestep > 0)
	    {
	      std::cout << "Redo time step. Accuracy not good enough.\n" << std::endl;
	    }

	  solution = tmp_solution_for_time_adaptivity;
	  time = tmp_time;
	  tmp_timestep = timestep;

	  // Large time step 2k
	  timestep = 2.0 * timestep;
	  old_timestep_solution = tmp_solution_for_time_adaptivity;
	  solve_spatial_problem ();

	  // Evaluate goal functional
	  double gf_stress_x_2k = goal_functional_stress_x ();



	  std::cout << "-------------------------" << std::endl;
	  std::cout << "Two small time steps\n" << std::endl;
	  // Two times small time steps
	  time = tmp_time;
	  timestep = tmp_timestep;
	  solution = tmp_solution_for_time_adaptivity;

	  for (unsigned int l=0;l<2;l++)
	    {
	      old_timestep_solution = solution;
	      solve_spatial_problem ();

	       time += timestep; 
	    }

	  // Evaluate goal functional
	  double gf_stress_x_1k = goal_functional_stress_x ();


	  // Error evaluation
	  // ABS 
	  double ABS_functional_error = 0.0;

	  ABS_functional_error= std::abs(gf_stress_x_2k - gf_stress_x_1k);
	  REL_functional_error = ABS_functional_error / std::abs(gf_stress_x_1k);

	  double gamma_safety_factor = 1.0;
	  double power_time_step_control = 0.06666666666666666666;
	  double theta_time_step_control = gamma_safety_factor * 
	    std::pow((TOL_for_time_adaptivity) / REL_functional_error, power_time_step_control); 
	  
	  if (use_time_step_control)
	    {
	      if (theta_time_step_control < 1.0 ||
		  1.2 < theta_time_step_control)
		{
		  timestep = theta_time_step_control * timestep;
		}
	      else 
		timestep = 1.0 * timestep;
	    }
	  
	// PI-version (John and Rang) - does not yield satisfactory results
		//double rho_time_adapt = 0.9;
		//timestep = rho_time_adapt * timestep * timestep / old_timestep * 
		//  std::sqrt(TOL_for_time_adaptivity * REL_old_timestep_functional_error / (REL_functional_error * REL_functional_error));
		
		// Growth, maximum and minimum timesteps
		if (timestep > (timestep_growth_factor * tmp_timestep))
		  timestep = timestep_growth_factor * tmp_timestep;
		
		if (timestep >= timestep_upper_bound) 
		  timestep = timestep_upper_bound;
		
		if (timestep <= timestep_lower_bound)
		  timestep = timestep_lower_bound;
		
		std::cout << std::endl;
		std::cout << "-------------------------" << std::endl;
		std::cout << "Functional:   " << std::setprecision(5) << time << "   " <<  std::setprecision(10) << gf_stress_x_1k << std::endl;

		std::cout << "TOL:          " << std::setprecision(5) << time << "   " << std::setprecision(10) << TOL_for_time_adaptivity << std::endl;
		std::cout << "ABS error:    " << std::setprecision(5) << time << "   " << std::setprecision(10) << ABS_functional_error << std::endl;
		std::cout << "REL error:    " << std::setprecision(5) << time << "   " << std::setprecision(10) << REL_functional_error << std::endl;
		std::cout << "New timestep: " << std::setprecision(5) << time << "   " << std::setprecision(10) << timestep << std::endl;
		std::cout << "Old timestep: " << std::setprecision(5) << time << "   " << std::setprecision(10) << tmp_timestep << std::endl;
		std::cout << std::endl;
		
		counter_redo_timestep += 1;
		
	      } while (timestep < timestep_rejection_factor * tmp_timestep);
	      
	REL_old_timestep_functional_error = REL_functional_error;
	old_timestep = tmp_timestep;


      }


      // Inside adaptive time loop now
      //time += timestep;
      //timestep = tmp_timestep;
      
	
      // Compute functional values
      // Evaluate the summed goal functionals with the 'old' timestep,
      // thus we need to copy from tmp_timestep
      double new_timestep = timestep;


      // Using time step control, we always advance three time steps.
      // Thus we need to multiply with '2' for the goal functional evaluations.
      // TODO: introduce variable for time_adaptivity
      if (time_stepping_scheme == "BE_adaptive" && timestep_number > 0)
	timestep = 2.0 * tmp_timestep;
      else 
	timestep = 1.0 * tmp_timestep;



      // Compute functional values
      std::cout << std::endl;
      if (test_case == "miehe_shear" ||
	  test_case == "miehe_tension" ||
	  test_case == "l_shaped")
	compute_functional_values();

      if (test_case == "Sneddon")
	compute_functional_values_Sneddon();

      if (test_case == "screw_domi")
	{
	  compute_stress_per_cell();
	  compute_energy();
	}

      if (test_case == "Sneddon3D")
	compute_cod_Sneddon3D();


      // Reset to standard timestep
      timestep = new_timestep;
      
      // Write solutions 
      // Plot corrector solution on new mesh (i.e., the final solution)
      if ((timestep_number % output_skip == 0))
	output_results (timestep_number+1,solution);
      

      if (bool_use_strain_history)
	{
	  // Update stress history
	  update_quadrature_point_history ();
	}

      //      if (timestep_number  == refine_mesh_1st ||
      //  timestep_number  == refine_mesh_2nd ||
      //  timestep_number  == refine_mesh_3rd ||
      //  timestep_number  == refine_mesh_4th ||
      //  timestep_number  == refine_mesh_5th
      //  )			      			      			     
	{
	  std::cout << "Refinement cycle " 
		    << refinement_cycle 
		    << "\n================== "
		    << std::endl;
	  
	  refine_mesh ();
	  ++refinement_cycle;		
	}


 
      

      
      ++timestep_number;

    }
  while ((timestep_number <= max_no_timesteps) && (time <= end_time_value));
  
  
}

// The main function looks almost the same
// as in all other deal.II tuturial steps. 
int main () 
{
  try
    {
      deallog.depth_console (0);

      const unsigned int dimension = 2;
      Dynamic_Fracture_Problem<dimension> flow_problem(1);
      flow_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}




