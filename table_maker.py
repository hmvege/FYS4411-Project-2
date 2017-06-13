import os, sys

blocking_output_directory = "blocking_output/"
terminal_output_directory = "terminal_output/"

def retrieve_file_contents(files):
	return_list = []
	for file in files: # Loops over all files in folder
		with open(file,"r") as f:
			for line in f: # Loops over all lines in file
				return_dict = {}
				vals = line.split()
				for i in range(len(vals)/2): # Loops over all values
					return_dict[vals[2*i]] = vals[2*i+1] # Puts values pair-wise into a dictionary
				return_list.append(return_dict)
	return return_list

def get_files(folders):
	return [[blocking_output_directory+folder+"/"+file for file in os.listdir(blocking_output_directory + folder) if file != ".DS_Store"] for folder in folders] # Returns all files in folder

def create_table(value_matrix, column_values=[], row_values=[]):
	print len(value_matrix[0])
	N_rows = len(value_matrix)
	if len(row_values) != 0: 
		N_rows += 1
		value_matrix.insert(0,row_values)
	N_cols = len(value_matrix[0])
	# if len(column_values) != 0: N_cols += 1
	tab = ""
	for i in xrange(N_rows):
		row = ""
		for j in xrange(N_cols):
			row += value_matrix[i][j]
			if j < N_cols-1: row += " &"
		row += " \\\ "
		tab += row + "\n"
	return tab

# Hard-coded
hc_folders = ["2e_jastrowWithCoulomb","2e_jastrow","2e_plain"]
hc_files = [f[0] for f in get_files(hc_folders)]
hc_data = retrieve_file_contents(hc_files)
print hc_data
var_params_header_string = []
energy_header_string_hc = [	r"$\langle E \rangle$", r"$\sigma^2$", r"", r"$\alpha$", r"$\beta$", r"Acceptance"]

hc_value_matrix = [[item for sublist in [["%15.8g" % float(energies["E"]),"%15.8g" % float(energies["EVar"])] for energies in hc_data] for item in sublist]]
print hc_value_matrix
print create_table(hc_value_matrix, row_values = energy_header_string_hc)

# No interaction
no_int_folders = ["no_interaction"]
no_int_files = [f[0] for f in get_files(no_int_folders)]
no_int_data = retrieve_file_contents(no_int_files)
# print no_int_data
energy_header_string_no_int = [r"$N$", r"$\omega$", 
							r"$\langle E_\text{J} \rangle$", r"$\sigma^2_{J}$",
							r"$\langle E_\text{J, no Coulomb} \rangle$", r"$\sigma^2_{J, no Coulomb}$",
							r"$\langle E_\text{no Coulomb} \rangle$", r"$\sigma^2_{no Coulomb}$",
							r"$Acceptance$"]


# No importance sampling
no_imp_folders = ["no_imp"]


# Importance sampling
imp_folders = ["imp"]
