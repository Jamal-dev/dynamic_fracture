#include "./../../include/dynamic_fracture.h"

// In this function, we solve the linear systems
// inside the nonlinear Newton iteration. We only
// use a direct solver from UMFPACK.
template <int dim>
void 
Dynamic_Fracture_Problem<dim>::solve () 
{
  #if DEAL_II_VERSION_GTE(9,0,0)
    timer.enter_subsection("Solve linear system.");
  #else
    timer.enter_section("Solve linear system.");
  #endif
  
  Vector<double> sol, rhs;    
  sol = newton_update;    
  rhs = system_rhs;
  
  SparseDirectUMFPACK A_direct;
  A_direct.factorize(system_matrix);     
  A_direct.vmult(sol,rhs); 
  newton_update = sol;
  
  // TODO: constraints.distrubte is it making sense here?
  constraints.distribute (newton_update);
  #if DEAL_II_VERSION_GTE(9,0,0)
    timer.leave_subsection("Solve linear system.");
  #else
    timer.exit_section();
  #endif
  
}
