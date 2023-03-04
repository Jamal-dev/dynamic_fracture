file_path = "example1/mesh_1/parameters_current_simulation.csv"
temp = file_path.split("/")
output_file_path = f'{temp[0]}/{temp[1]}/' + 'latex_parameters_setting.txt'
results = {}
with open(file_path,'r') as f:
	content = f.readlines ()
	for row in content:
		values = row.split('\n')[0].split(',')[:-1]
		for key, value in zip(values[0::2],values[1::2]):
			results[key] = value
max_no_steps = 1.3e4
with open(output_file_path,'w') as f:
	f.writelines(f''' \t Cells & {results['Cells']} & DoFs & {results['DOFs']} & ~ & ~ \\\\ \hline
    DoFs\_u & {results['Dofs_u']} & DoFs\_p & {results['Dofs_p']} & DoFs\_v & {results['Dofs_v']} \\\\ \hline
    $\lambda$ & \\numS{{{results['lambda']}}} \si{{\pascal}} & $\mu$ & \\numS{{{results['mu']}}} \si{{\pascal}} & ~ & ~\\\\ \hline
    $\\rho$ & \\num{{{results['density']}}} & $G_c$ & \\num{{{results['G_c']}}} & min\_dia & {float(results['eps'])/2} \\\\ \hline
    $E$ & \\numS{{{results['E']}}} \si{{\pascal}} & $\\nu$ & {results['nu']} & $l_0$ & {results['eps']}\\\\ \hline
    Time steps & \\num{{{results['timestep']}}} & Final time & \\num{{{max_no_steps*float(results['timestep'])}}} & ~ & ~\\\\ \hline ''')
print(f'File has written successfully... at')
print(output_file_path)
