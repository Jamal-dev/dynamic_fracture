#!/usr/bin/env python
# coding: utf-8

# In[12]:


import re
from path import Path





def write_code(modify=False, path=None, match_str = None, target_str = None, add_str = None):
    # match_str = f'else if (current_test_case == test_cases::{source.upper()})'
    # target_str = f'else if (current_test_case == test_cases::{target.upper()})'
    source_file = Path(path)
    if not modify:
        with open(path,'r') as f:
            lines = f.readlines()
            for i,line in enumerate(lines):
                if target_str in line:
                    print(f'file: {path}, has already function: {target_str}')
                    return 
    with open(path, 'r') as f:
        lines = f.readlines()
        flag = False
        copy_lines = []
        count_left_bracket = 0
        count_right_bracket = 0
        for i, line in enumerate(lines):
            if match_str in line:
                flag = True
            if flag :
                if "}" in line:
                    count_right_bracket += 1
                    if count_left_bracket == count_right_bracket:
                        flag = False
                        copy_lines.append("\n\t}\n")
                        break
                    else:
                        copy_lines.append(line)
                        continue
                    
                else:
                    if '{' in line:
                        count_left_bracket += 1
                        if count_left_bracket == 1:
                            copy_lines.append('\t'+line)
                        else:
                            copy_lines.append(line)
                    else:
                        copy_lines.append(line)

        
        # Insert the new line
        lines.insert(i, f'\t}}\n{target_str}\n\t{{\n')
        starting_index = i + 1
        output_file = path
            
        for j,line in enumerate(copy_lines[2:-1]):
            lines.insert(i+j+1, line)

        
    print(f'Function: {target_str} added to file: {output_file}')
    print(f'Edit this function at line: {starting_index}')
    with open(path, 'w') as f:
        f.writelines(lines)
        

    
    







# # Write test cases

# In[62]:


def write_test_cases(modify=False, path=None, match_str = None, target_name = None, add_str = None,tabs=1,is_cond=False):
    if not modify:
        with open(path,'r') as f:
            lines = f.readlines()
            flag = False
            for i,line in enumerate(lines):
                if match_str in line:
                    flag = True
                if flag and target_name in line:
                    print(f'file: {path}, has already function: {target_name}')
                    return 
    with open(path, 'r') as f:
        lines = f.readlines()
        flag = False
        for i, line in enumerate(lines):
            if match_str in line:
                flag = True
            if flag :
                if "}" in line:
                    flag = False
                    tabs_ = ''
                    for _ in range(tabs):
                        tabs_ += '\t'
                    if is_cond:
                        name_field = f'{tabs_}constexpr bool Is_{target_name.lower()}() const {{ return value == {target_name.upper()}; }}'
                        lines.insert(i+1, name_field)
                    else:
                        name_field = f'{tabs_}, {target_name.upper()}\n'
                        lines.insert(i, name_field)
                    break
                

        
        
        # Insert the new line
        
        starting_index = i + 1
        output_file = path
        
        
    print(f'Function: {target_name} added to file: {output_file}')
    print(f'Edit this function at line: {starting_index}')
    with open(path, 'w') as f:
        f.writelines(lines)
        if add_str is not None:
            f.write(add_str)



# In[60]:


# goal functional
def write_goal_functional(modify=False, path=None, match_str = None, target_name = None, add_str = None,tabs=1):
    if not modify:
        with open(path,'r') as f:
            lines = f.readlines()
            for i,line in enumerate(lines):
                if target_name in line:
                    print(f'file: {path}, has already function: {target_name}')
                    return 
    with open(path, 'r') as f:
        lines = f.readlines()
        flag = False
        for i, line in enumerate(lines):
            if match_str in line:
                flag = True
            if flag :
                tabs_ = ''
                for _ in range(tabs):
                    tabs_ += '\t'
                name_field = f'{tabs_} {target_name} ||\n'
                lines.insert(i+1, name_field)
                print(f'Added line: {i+2}')
                break
                

        
        
        # Insert the new line
        
        starting_index = i + 1
        output_file = path
        
        
    print(f'Function: {target_name} added to file: {output_file}')
    print(f'Edit this function at line: {starting_index}')
    with open(path, 'w') as f:
        f.writelines(lines)
        if add_str is not None:
            f.write(add_str)



# In[51]:


def write_parameters(previous_file_modify=False, path=None, parameters:dict = None,mesh_path=''):
    source_file = Path(path)
    target_file = source_file.parent / Path(f'{target.lower()}.cpp')
    source_name_function = f'void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_{source.lower()} ()'
    print(f'Adding changes to file: {target_file}')
    target_name_function = f'void Dynamic_Fracture_Problem<dim>::set_runtime_parameters_{target.lower()} ()'
    target_name_function_only = f'set_runtime_parameters_{target.lower()} ();'
    def check_if_already_changed(path,target_line):
        with open(path,'r') as f:
            lines = f.readlines()
            for i,line in enumerate(lines):
                if target_line in line:
                    return False
        return True
    if not previous_file_modify:
        if target_file.exists():
            print(f'file: {target_file}, already exists')
            return
         
    with open(path, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if source_name_function in line:
                lines[i] = target_name_function
                print(f'Line: {i+1} modified')
            if 'lame_coefficient_mu' in  line or 'lame_coefficient_lambda' in line:
                break
            for key in parameters.keys():
                if key in line:
                    lines[i] = f'\t{key} = {parameters[key]};\n'
                    print(f'Line: {i} modified')

    with open(target_file, 'w') as f:
        f.writelines(lines)
    with open(target_file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 'std::string grid_name;' in line:
                continue
            if "grid_name" in line:
                lines[i] = f'\t\tgrid_name = "{mesh_path}";\n'
                print(f'Line: {i} modified')
                break
    with open(target_file, 'w') as f:
        f.writelines(lines)
    with open(target_file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "std::string parent_dir" in line:
                lines[i] = f'\tstd::string parent_dir = "./results/patrick_tests/{target.lower()}";\n'
                print(f'Line: {i} modified')
                break
    with open(target_file, 'w') as f:
        f.writelines(lines)
    with open(target_file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 'filename_basis  = parent_dir + "/"' in line:
                lines[i] = f'\tfilename_basis  = parent_dir + "/" +  "solution_{target.lower()}_test_";\n'
                print(f'Line: {i} modified')
                break
                    
          
    with open(target_file, 'w') as f:
        f.writelines(lines)

    target_line = f'#include "source/parameters/{target.lower()}.cpp"'
    target_file = 'dynamic_fracture.cc'
    if check_if_already_changed(target_file,target_line):
        print('Adding parameters to {target_file} file')
        with open(target_file,'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if '#include "source/parameters/p_mesh_1.cpp"' in line:
                    lines.insert(i+1, f'{target_line}\n')
                    print(f'Line: {i+2} Added')
                    break
        with open(target_file,'w') as f:
            f.writelines(lines)
    target_line = f'else if (current_test_case == test_cases::{target.upper()})'
    target_line_add = f'    set_runtime_parameters_{target}();'
    if check_if_already_changed(target_file,target_line):
        print(f'Adding parameters to {target_file} file')
        with open(target_file,'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if 'else if (current_test_case == test_cases::P_MESH_1)' in line:
                    lines.insert(i+2, f'\t{target_line}\n')
                    lines.insert(i+3, f'{target_line_add}\n')
                    print(f'Line: {i+3} Added')
                    print(f'Line: {i+4} Added')
                    break
        with open(target_file,'w') as f:
            f.writelines(lines)
    target_line = f'void {target_name_function_only}'
    target_file = 'include/dynamic_fracture.h'
    if check_if_already_changed(target_file,target_line):
        print(f'Adding parameters to {target_file} file')
        with open(target_file,'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if 'void set_runtime_parameters_p_mesh_1 ();' in line:
                    lines.insert(i+1, f'\tvoid {target_name_function_only}\n')
                    print(f'Line: {i+2} Added')
                    break
        with open(target_file,'w') as f:
            f.writelines(lines)
    
        



# In[52]:


def section(text):
    t = ''
    for i in range(len(text)+4):
        t += '='
    print(t)
    print(f'== {text} ==')
    print(t)
def calling_functions(parameters,mesh_path):
    section('Adding parameter file')
    path = f'source/parameters/{source.lower()}.cpp'
    write_parameters(previous_file_modify=True, 
                        path=path, 
                        parameters = parameters,
                        mesh_path = mesh_path)
    section('Adding test case into goal functional')
    path = 'source/compute_quantities/goal_functional_stress.cpp'
    match_str = 'if (current_test_case.Is_miehe_shear()  ||'
    target_name = f'current_test_case.Is_{target.lower()}()'
    write_goal_functional(modify = False,
                    path = path,
                    match_str = match_str,
                    target_name = target_name,
                    tabs = 3)
    section('Adding test case into dynamic fracture')
    path = 'include/test_cases.h'
    match_str = 'enum Value : uint8_t'
    target_name = target.upper()
    write_test_cases(modify = False,
                    path = path,  
                    match_str = match_str,
                    target_name = target_name)
    match_str = 'std::set<test_cases> all_cases()  {'
    write_test_cases(modify = False,
                    path = path,  
                    match_str = match_str,
                    target_name = target_name,
                    tabs = 5)
    match_str = 'Is_pmesh1()'
    write_test_cases(modify = False,
                    path = path,  
                    match_str = match_str,
                    target_name = target_name,
                    tabs = 2,
                    is_cond = True)
    section("Newton BC's")
    path = 'source/boundary_conditions/newton_bc.cpp'
    match_str = f'else if (current_test_case == test_cases::{source.upper()})'
    target_str = f'else if (current_test_case == test_cases::{target.upper()})'
    write_code(modify = False,
            path = path,
            match_str = match_str,
            target_str = target_str)
    section("Initial guess newton")
    # initial guess newton
    path = 'source/boundary_conditions/initial_guess_newton.cpp'
    match_str = f'else if (current_test_case == test_cases::{source.upper()})'
    target_str = f'else if (current_test_case == test_cases::{target.upper()})'
    write_code(modify = True,
            path = path,
            match_str = match_str,
            target_str = target_str) 
    section("Nonohomogenous direchlet BC's values")
    # nonohomogenu dirchlet boundary values
    path = 'source/boundary_conditions/nonhomoDirichletBoundary.cpp'
    match_str = f'else if (_test_case == test_cases::{source.upper()})'
    target_str = f'\telse if (_test_case == test_cases::{target.upper()})'
    write_code(modify = True,
            path = path,
            match_str = match_str,
            target_str = target_str,
            )
    section("Nonohomogenous direchlet BC's velocity values")
    # nonohomogenu dirchlet boundary velicity values
    path = 'source/boundary_conditions/nonhomoDirichletBoundaryVelocities.cpp'
    match_str = f'else if (_test_case == test_cases::{source.upper()})'
    target_str = f'else if (_test_case == test_cases::{target.upper()})'
    write_code(modify = True,
            path = path,
            match_str = match_str,
            target_str = target_str,
            )

    

def main():
    global target, source
    target = 'p_notched_cavity'
    source = 'p_mesh_1'
    parameters = dict()
    parameters['G_c'] = 8e3
    parameters['density_structure'] = 1.2e3
    parameters['poisson_ratio_nu'] = 1/3
    parameters['E_modulus'] = 72000e6
    parameters['timestep'] = 1e-7
    parameters['end_time_value'] = 1e-1
    mesh_path = f'mesh_files/example2/{target.lower()}.msh'
    calling_functions(parameters,mesh_path)

if __name__ == '__main__':
    main()