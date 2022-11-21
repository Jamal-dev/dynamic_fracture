First of all write your case in file source/test_cases.cppp
Create your file in source/parameters/set_runtime_parameters_yourcase();
Add this function inside
include/dynamic_fracture.h

Then add your test case in 
source/boundary_conditions/initial_bc.cpp
source/boundary_conditions/newton_bc.cpp
source/boundary_conditions/nonhomoDirichletBoundary.cpp

Inside source/dynamic_fracture.cc
Go to line 
"// Setting specific parameters for each problem"
Add your test case and aslo your 
set_runtime_parameters_yourcase();

Add goal_functional.cpp
if (current_test_case.IsMiehe_Shear()  || 
       current_test_case.IsMiehe_Tension()
       || current_test_case.IsPMesh1())