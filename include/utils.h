#ifndef utils_local_h
#define utils_local_h
#include <vector>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
using namespace dealii;





// First, we define tensors for the solution variables
// Note: note all of these definitions are needed in 
// phase-field fracture. Some of them come from 
// fluid-structure interaction computations.
namespace ALE_Transformations
{    

// get_grad_pf
 template <int dim> 
   inline
   Tensor<1,dim> 
   get_grad_pf (unsigned int q,
	       std::vector<std::vector<Tensor<1,dim> > > old_solution_grads)	 
   {     
     Tensor<1,dim> grad_pf;     
     grad_pf[0] =  old_solution_grads[q][dim][0];
     grad_pf[1] =  old_solution_grads[q][dim][1];
     if (dim == 3)
       grad_pf[2] = old_solution_grads[q][dim][2];
      
     return grad_pf;
   }

 
 // get_grad_u
 template <int dim> 
   inline
   Tensor<2,dim> 
   get_grad_u (unsigned int q,
	       std::vector<std::vector<Tensor<1,dim> > > old_solution_grads)	 
   {   
     Tensor<2,dim> grad_u;
    grad_u[0][0] =  old_solution_grads[q][0][0];
    grad_u[0][1] =  old_solution_grads[q][0][1];

    grad_u[1][0] =  old_solution_grads[q][1][0];
    grad_u[1][1] =  old_solution_grads[q][1][1];
    if (dim == 3)
      {
        grad_u[0][2] =  old_solution_grads[q][0][2];

        grad_u[1][2] =  old_solution_grads[q][1][2];

        grad_u[2][0] =  old_solution_grads[q][2][0];
        grad_u[2][1] =  old_solution_grads[q][2][1];
        grad_u[2][2] =  old_solution_grads[q][2][2];
      }

    return grad_u;
   }

// get_Identity
  template <int dim> 
    inline
    Tensor<2,dim> 
    get_Identity ()
    {   
      Tensor<2,dim> identity;
      identity[0][0] = 1.0;
      identity[1][1] = 1.0;
      if (dim == 3)
	identity[2][2] = 1.0;
            
      return identity; 
    }


// get_u
template <int dim> 
 inline
 Tensor<1,dim> 
 get_u (unsigned int q,
	std::vector<Vector<double> > old_solution_values)
   {
     Tensor<1,dim> u;     
     u[0] = old_solution_values[q](0);
     u[1] = old_solution_values[q](1);
     if (dim == 3)
       u[2] = old_solution_values[q](2);
     
     return u;          
   }

// get_u_LinU
 template <int dim> 
   inline
   Tensor<1,dim> 
   get_u_LinU (const Tensor<1,dim> phi_i_u)
   {
     Tensor<1,dim> tmp;     
     tmp[0] = phi_i_u[0];
     tmp[1] = phi_i_u[1];
     if (dim == 3)
       tmp[2] = phi_i_u[2];
     
     return tmp;    
   }

// get_v
 template <int dim> 
 inline
 Tensor<1,dim> 
 get_v (unsigned int q,
	std::vector<Vector<double> > old_solution_values)
   {
     Tensor<1,dim> v;     
     v[0] = old_solution_values[q](dim+1);
     v[1] = old_solution_values[q](dim+2);
     if (dim == 3)
       v[2] = old_solution_values[q](dim+3);
     
     return v;          
   }
 
// get_grad_v
 template <int dim> 
   inline
   Tensor<2,dim> 
   get_grad_v (unsigned int q,
	       std::vector<std::vector<Tensor<1,dim> > > old_solution_grads)	 
   {   
     Tensor<2,dim> grad_u;
    grad_u[0][0] =  old_solution_grads[q][dim+1][0];
    grad_u[0][1] =  old_solution_grads[q][dim+1][1];

    grad_u[1][0] =  old_solution_grads[q][dim+2][0];
    grad_u[1][1] =  old_solution_grads[q][dim+2][1];
    if (dim == 3)
      {
        grad_u[0][2] =  old_solution_grads[q][dim+1][2];

        grad_u[1][2] =  old_solution_grads[q][dim+2][2];

        grad_u[2][0] =  old_solution_grads[q][dim+3][0];
        grad_u[2][1] =  old_solution_grads[q][dim+3][1];
        grad_u[2][2] =  old_solution_grads[q][dim+3][2];
      }

    return grad_u;
   }


// get_deviator_norm
  template <int dim> 
  inline
  double
  get_deviator_norm (const Tensor<2,dim> deviator)	 
  {    
    if (dim == 2)
      {
	return std::sqrt(deviator[0][0] * deviator[0][0] 
			 + deviator[0][1] * deviator[0][1]
			 + deviator[1][0] * deviator[1][0]
			 + deviator[1][1] * deviator[1][1]);
      }
    else if (dim == 3)
      {
	// TODO: needs to be done
	return 0.0;


      }
    
  }
  
 // get diverence of u
 template <int dim>
  inline
  double
  get_divergence_u (const Tensor<dim,dim> grad_u)
  {

    if (dim == 2)
      {
        return grad_u[0][0] + grad_u[1][1];
      }
    else if (dim == 3)
      {
        return grad_u[0][0] + grad_u[1][1] + grad_u[2][2];
      }
     else
         std::logic_error("Function not yet implemented");


  }


}




// structure terms in ALE
// In the third namespace, we summarize the 
// constitutive relations for the structure equations.
namespace Structure_Terms_in_ALE
{
  // Green-Lagrange strain tensor
  template <int dim> 
  inline
  Tensor<2,dim> 
  get_E (const Tensor<2,dim> F_T,
	 const Tensor<2,dim> F,
	 const Tensor<2,dim> Identity)
  {    
    return 0.5 * (F_T * F - Identity);
  }

  template <int dim> 
  inline
  double
  get_tr_E (const Tensor<2,dim> E)
  {     
    return trace (E);
  }

  template <int dim> 
  inline
  double
  get_tr_E_LinU (unsigned int q, 
		 const std::vector<std::vector<Tensor<1,dim> > > old_solution_grads,
		 const Tensor<2,dim> phi_i_grads_u)	    
  {
    return ((1 + old_solution_grads[q][0][0]) *
	    phi_i_grads_u[0][0] + 
	    old_solution_grads[q][0][1] *
	    phi_i_grads_u[0][1] +
	    (1 + old_solution_grads[q][1][1]) *
	    phi_i_grads_u[1][1] + 
	    old_solution_grads[q][1][0] *
	    phi_i_grads_u[1][0]); 
  }
  
}



// Bitmap class
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

/*
    Bitmap function class
*/

template <int dim>
class BitmapFunction : public Function<dim>
{
public:
  BitmapFunction(const std::string &filename,
                 double x1_, 
                 double x2_, 
                 double y1_, 
                 double y2_, 
                 double minvalue_, 
                 double maxvalue_);

  virtual double value (const Point<dim> &p,
                const unsigned int /*component*/ = 0) const;
private:
  BitmapFile f;
  double x1,x2,y1,y2;
  double minvalue, maxvalue;
};
///////////////////////////////////////////////////////////////



#endif
