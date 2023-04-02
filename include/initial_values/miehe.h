#ifndef miehe_initial_values_local_h
#define miehe_initial_values_local_h
#include <deal.II/lac/vector.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include "./../test_cases.h"
#include <fstream>
#include <math.h>
using namespace dealii;
using namespace std;
struct PointStruct{
  double x;
  double y;
};
struct LineStruct{
  // model a x + by + c = 0
  double a;
  double b;
  double c;
};

template <int dim>
  class InitialValuesPhaseField : public Function<dim>
  {
    public:
      InitialValuesPhaseField (const double & alpha_eps, 
			  const bool & bool_initial_crack_via_phase_field, const test_cases & test_case) : Function<dim>(dim+1+dim) 
    {
      _alpha_eps = alpha_eps;
      _bool_initial_crack_via_phase_field = bool_initial_crack_via_phase_field;
      _test_case = test_case;
      _min_radius_mesh = _alpha_eps/2.0; // assuming diameter of the mesh; although it's 2 or 3 times the diameter
      _min_radius_mesh_squared = _min_radius_mesh*_min_radius_mesh;
      width = 10;
      height = 10;
      eps_mesh = 1e-6;
      PointStruct mid_point = {width/2,height/2};
      PointStruct p1 = {0,mid_point.y-eps_mesh};
      PointStruct p2 = {0,mid_point.y+eps_mesh};
      m1 = (mid_point.y - p1.y)/(mid_point.x - p1.x);
      m2 = (mid_point.y - p2.y)/(mid_point.x - p2.x);
      b1 = mid_point.y - m1*mid_point.x;
      b2 = mid_point.y - m2*mid_point.x;
      notched_tip = {width/2,height/2};
      before_after = false;

      if (_test_case.Is_p_notched_cavity())
        {
          width = 10;
          height = 10;
          double notched_vert_dist_from_top_edge = 1;
          double notched_width = 1.5;
          double notched_height = 3;

          PointStruct notched_upper = {width, height- notched_vert_dist_from_top_edge};
          PointStruct notched_mid = {notched_upper.x -  notched_width, notched_upper.y -  notched_height/2.0};
          PointStruct notched_lower = {notched_upper.x, notched_mid.y -  notched_height/2.0};

          notched_tip = notched_mid;

          m1 = (notched_mid.y - notched_upper.y)/(notched_mid.x - notched_upper.x);
          m2 = (notched_mid.y - notched_lower.y)/(notched_mid.x - notched_lower.x);
          b1 = notched_mid.y - m1 * notched_mid.x;
          b2 = notched_mid.y - m2 * notched_mid.x;
          before_after = true;
          
          
        }
        // a x + b y + c = 0
        // dy/dx = -a/b = m1  -> a = -b*m1
        // b1 = -c/b  => 
        // a*a + b*b = 1
        // b = 1/ sqrt(1+m1*m1)
        // a = -b*m1
        // c = -b1*b 
        line1.b = 1/sqrt(m1*m1 + 1);
        line1.a = -line1.b*m1;
        line1.c = -b1*line1.b;
        line2.b = 1/sqrt(m2*m2 + 1);
        line2.a = -line2.b*m2;
        line2.c = -b2*line2.b;
      
    }

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;
      bool check_near_line(const double x, const double y,  const LineStruct &line) const;
      bool check_near_tip(const double x, const double y) const;
      bool region_before_notch(const double x, const double y) const;
      bool region_after_notch(const double x, const double y) const;
      bool select_region(const double x, const double y, const bool before_after) const;

  private:
    double _alpha_eps;
    bool _bool_initial_crack_via_phase_field;
    test_cases _test_case;
    bool before_after;
    double width;
    double height;
    double eps_mesh;
    double m1,m2,b1,b2;
    double _min_radius_mesh;
    double _min_radius_mesh_squared;
    LineStruct line1;
    LineStruct line2;
    PointStruct notched_tip;

  };

  template <int dim>
  bool
  InitialValuesPhaseField<dim>::check_near_tip(const double x, const double y) const
  {
    double dist = (x-notched_tip.x)*(x-notched_tip.x) + (y-notched_tip.y)*(y-notched_tip.y);
    if (dist < _min_radius_mesh_squared)
      return true;
    else
      return false;
  }

  template <int dim>
  bool
  InitialValuesPhaseField<dim>::region_before_notch(const double x, const double y) const
  {
    if (x< notched_tip.x-_min_radius_mesh || y< notched_tip.y-_min_radius_mesh || y> notched_tip.y+_min_radius_mesh)
      return true;
    else
      return false;
  }

  template <int dim>
  bool
  InitialValuesPhaseField<dim>::region_after_notch(const double x, const double y) const
  {
    if (x> notched_tip.x+_min_radius_mesh || y< notched_tip.y-_min_radius_mesh || y> notched_tip.y+_min_radius_mesh)
      return true;
    else
      return false;
  }

  template <int dim>
  bool
  InitialValuesPhaseField<dim>::select_region(const double x, const double y, const bool before_after) const
  {
    if (before_after)
      return region_before_notch(x,y);
    else
      return region_after_notch(x,y);
  }

  template <int dim>
  bool
  InitialValuesPhaseField<dim>::check_near_line(const double x, const double y, const LineStruct &line) const
  {
    if (select_region(x,y,before_after))
      return false;    
    double dist = fabs(line.a*x + line.b*y + line.c)/sqrt(line.a*line.a + line.b*line.b);
    if (dist < _min_radius_mesh)
      return true;
    else
      return false;
  }


  template <int dim>
  double
  InitialValuesPhaseField<dim>::value (const Point<dim>  &p,
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
          cout<<"Inside of _bool_initial_crack_via_phase_field"<<endl;
          return 0.0; 

        }
      else 
        return 1.0;
          }
      }
	else // _bool_initial_crack_via_phase_field == false
      {
        if (component == 2)   
          {
            // fstream myfile;
            // myfile.open ("testing_coord.txt", ios::out | ios::app);
            // myfile<<p(0)<<endl;
            // myfile<<p(1)<<endl;
            // myfile.close();
            if ( _test_case.Is_p_notched_cavity())
              {
                
                // checking if the point lies at the crack line
                //if (p(1) == m1*p(0) + b1 || p(1) == m2*p(0) + b2)
                //if (check_near_line(p(0),p(1),line1) || check_near_line(p(0),p(1),line2))
                if (check_near_tip(p(0),p(1)))
                  {
                    // cout<<"******************************************************"<<endl;
                    // cout<<"p(0) = "<<p(0)<<endl;
                    // cout<<"p(1) = "<<p(1)<<endl;
                    // cout<<"m1*p(0) + b1 ="<< m1*p(0) + b1<<endl;
                    // cout<<"m2*p(0) + b2 ="<< m2*p(0) + b2<<endl;
                    // cout<<"******************************************************"<<endl;
                    return 0.0;
                    }
                else
                  {
                    // cout<<"Not on the notched lines"<<endl;
                    return 1.0;}

              }            
            else
              return 1.0;
          }

      }

    // for all other components are 0
    return 0.0;


  }


  template <int dim>
  void
  InitialValuesPhaseField<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int comp=0; comp<this->n_components; ++comp)
      values (comp) = InitialValuesPhaseField<dim>::value (p, comp);
  }

#endif // miehe_initial_values_local_h