import os, sys, re

# Global constants
blocking_output_directory = "blocking_output/"
terminal_output_directory = "terminal_output/"
n_precision = 5
box_width = 2012122

def retrieve_file_contents(files,single_data_file=False):
	return_list = []
	for file in files: # Loops over all files in folder
		with open(file,"r") as f:
			if single_data_file: # If we have a single file, and all data is contained in that one
				return_dict = {}
			for line in f: # Loops over all lines in file
				if not single_data_file:
					return_dict = {}
				vals = line.split()
				for i in range(len(vals)/2): # Loops over all values
					return_dict[vals[2*i]] = vals[2*i+1] # Puts values pair-wise into a dictionary
				if not single_data_file:
					return_list.append(return_dict)
			if single_data_file:
				return_list.append(return_dict)
	return return_list

def get_files(folders, base):
	return [[base+folder+"/"+file for file in os.listdir(base + folder) if file != ".DS_Store"] for folder in folders] # Returns all files in folder

def create_table(value_matrix, column_values=[], row_values=[]):
	N_cols = len(value_matrix[0])
	if len(column_values) != 0:
		N_cols += 1
		for i in xrange(len(value_matrix)):
			value_matrix[i].insert(0,column_values[i])		
	N_rows = len(value_matrix)
	if len(row_values) != 0: 
		N_rows += 1
		value_matrix.insert(0,row_values)
	tab = ""
	print value_matrix
	for i in xrange(N_rows):
		row = ""
		for j in xrange(N_cols):
			row += value_matrix[i][j]
			if j < N_cols-1: row += " &"
		row += " \\\ "
		tab += row + "\n"
	return tab

clean_string = lambda s : "%15.5g" % float(s)
getDigit = lambda s : float(re.findall(r"[\d]+",s)[0])
getFraction = lambda s : float(re.findall(r"[\d.]+",s)[0])

def create_value_matrix(block_data, terminal_data,item_list):
	value_matrix = []
	for bd,td in zip(block_data,terminal_data):
		line = []
		# bd = block_data[i]
		# td = terminal_data[i]
		if "N" 		in item_list:			line.append(clean_string(bd["N"]))
		if "Omega" 	in item_list: 			line.append(clean_string(bd["Omega"]))
		if "E" 		in item_list:			line.append(clean_string(bd["E"]))
		if "EVar" 	in item_list:	 		line.append(clean_string(bd["EVar"]))
		if "KE" 	in item_list: 			line.append(clean_string(td["KE"]))
		if "PE" 	in item_list: 			line.append(clean_string(td["PE"]))
		if "Alpha" 	in item_list: 			line.append(clean_string(bd["Alpha"]))
		if "Beta" 	in item_list: 
			try: 
				line.append(clean_string(bd["Beta"]))
			except KeyError:
				line.append(" "*15)
		if "N_MC" 	in item_list: 			line.append("%15g" % int(bd["N_MC"]))
		if "AcceptanceRate" in item_list: 	line.append("%15.3g" % float(td["AcceptanceRate"]))
		value_matrix.append(line)
	return value_matrix

def particle_info_getter(line):
	line = line.split("/")[-1]
	line = line.split("_")
	particle_number = [item for item in line if "Particle" in item]
	particle_omega = [item for item in line if "omega" in item]
	return [getDigit(particle_number[0]),getFraction(particle_omega[0])]

#sort_after_omega
def sort_terminal_files(files):
	# Sorts after particle, then after omega
	new_list = []
	for file in files:
		new_list.append([file,particle_info_getter(file)])
	new_list.sort(key = lambda x: (x[1][0],-x[1][1]))
	for i in range(len(new_list)):
		files[i] = new_list[i][0]
	return files

# Cases ==========================================================================================================================================
# Hard-coded
hc_folders = ["2e_jastrowWithCoulomb","2e_jastrow","2e_plain"]
hc_blocking_files = [f[0] for f in get_files(hc_folders,blocking_output_directory)]
hc_terminal_files = [f[0] for f in get_files(hc_folders,terminal_output_directory)]
hc_blocking_data = retrieve_file_contents(hc_blocking_files)
hc_terminal_data = retrieve_file_contents(hc_terminal_files,True)

hc_header_string = [r"Case",r"$\langle E \rangle$", r"$\sigma^2$", r"$\langle K\rangle$", r"$\langle V \rangle$", r"$\alpha$", r"$\beta$", r"$N_{MC}$", r"Acceptance"]
hc_col_string = ["Jastrow","Jastrow, no Coulomb", "No Jastrow, no Coulomb"]

hc_value_matrix = []
for i in range(len(hc_folders)):
	# print hc_blocking_data[i]
	# print hc_terminal_data[i], "\n"
	line = []
	bd = hc_blocking_data[i]
	td = hc_terminal_data[i]
	line.append(clean_string(bd["E"]))
	line.append(clean_string(bd["EVar"]))
	line.append(clean_string(td["KE"]))
	line.append(clean_string(td["PE"]))
	line.append(clean_string(bd["Alpha"]))
	try:
		line.append(clean_string(bd["Beta"]))
	except KeyError:
		line.append(" "*15)
	line.append("%15g" % int(bd["N_MC"]))
	line.append("%15.3g" % float(td["AcceptanceRate"]))
	hc_value_matrix.append(line)

# print "Hard-coded table:\n"
# print create_table(hc_value_matrix, row_values=hc_header_string, column_values=hc_col_string), "\n\n"

# No interaction
no_int_folders = ["no_interaction"]
no_int_blocking_files = get_files(no_int_folders,blocking_output_directory)[0]
no_int_terminal_files = get_files(no_int_folders,terminal_output_directory)[0]
no_int_terminal_files = sort_terminal_files(no_int_terminal_files)
no_int_blocking_data = retrieve_file_contents(no_int_blocking_files)
no_int_terminal_data = retrieve_file_contents(no_int_terminal_files,True)
no_int_header_string = [r"$N$", 
						r"$\omega$", 
						r"$\langle E \rangle$", 
						r"$\sigma^2$", 
						r"$\langle K\rangle$", 
						r"$\langle V \rangle$", 
						r"$N_{MC}$",
						r"Acceptance"]
# item_to_include_list = ["N","Omega","E","EVar","KE","PE","Alpha","Beta","N_MC","AcceptanceRate"]

no_int_value_matrix = create_value_matrix(no_int_blocking_data,no_int_terminal_data,["N","Omega","E","EVar","KE","PE","N_MC","AcceptanceRate"])
# print no_int_value_matrix
print "No interaction table"
print create_table(no_int_value_matrix, row_values=no_int_header_string)

# No importance sampling
no_imp_folders = ["no_imp"]


# Importance sampling
imp_folders = ["imp"]
