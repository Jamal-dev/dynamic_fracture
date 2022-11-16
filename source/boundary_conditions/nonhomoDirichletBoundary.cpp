#include "./../../include/boundary_conditions/nonhomoDirichletBoundary.h"

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
  else if (_test_case == "dynamic_slit")
    {
      // Dynamic slit
	  // Important: dynamic slit is not from 0 to 1
	  // max(y) = 0.04, max(x) = 0.1
	  const float max_y = 0.04;
	  const float max_x = 0.1;
      if (component == comp.disp_y ) // u_y
		{
		return ( ((p(1) == max_y) && (p(0) <= max_x) && (p(0) >= 0.0)) 
			? 			
			(1.0) * _time * dis_step_per_timestep : 0 ); 
		
		}
    }
  else if (_test_case == "miehe_shear")
    {
      // Miehe shear
      if (component == 0) // u_x
		{
		return ( ((p(1) == 1.0) )
			?
			(-1.0) * _time *dis_step_per_timestep : 0 );
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
