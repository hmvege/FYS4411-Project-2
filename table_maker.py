import os, sys, re

# Global constants
blocking_output_directory = "blocking_output/"
terminal_output_directory = "terminal_output/"
n_precision = 5

#n 2,6,12,20
#omega 1,0.5,0.1,0.05,0.01
ERef = ["3.00000","1.659722",0,0,0,
		"20.1737", "11.8055",0,0,0,
		"65.7409", "39.2194",0,0,0,
		"155.9601","93.9891",0,0,0]

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
			value_matrix[i].insert(0,"%15s" % column_values[i])		
	N_rows = len(value_matrix)
	if len(row_values) != 0: 
		N_rows += 1
		value_matrix.insert(0,["%15s" % i for i in row_values])
	tab = ""
	for i in xrange(N_rows):
		row = ""
		for j in xrange(N_cols):
			if type(value_matrix[i][j]) == list:
				print value_matrix[i][j]
				raw_input(2)
			row += value_matrix[i][j]
			if j < N_cols-1: row += " &"
		row += " \\\ "
		if i==0: row += "\hline"
		tab += row + "\n"
	return tab

clean_string = lambda s : r"$%15.5g$" % float(s)
clean_omega = lambda s : r"$%15g$" % float(s)
getDigit = lambda s : float(re.findall(r"[\d]+",s)[0])
getFraction = lambda s : float(re.findall(r"[\d.]+",s)[0])

def clean_eref(s):
	if s == 0:
		return "-"+14*" " 
	else:
		return clean_string(s)

def create_value_matrix(block_data, terminal_data,item_list):
	value_matrix = []
	for bd,td in zip(block_data,terminal_data):
		line = []
		# bd = block_data[i]
		# td = terminal_data[i]
		if "N" 		in item_list:			line.append(clean_string(bd["N"]))
		if "Omega" 	in item_list: 			line.append(clean_omega(bd["Omega"]))
		if "E" 		in item_list:			line.append(clean_string(bd["E"]))
		if "EwoJ"	in item_list:			line.append(clean_string(bd["EwoJ"]))
		if "ERef"	in item_list:			line.append(clean_string(bd["ERef"]))
		if "EVar" 	in item_list:	 		line.append(clean_string(bd["EVar"]))
		if "KE" 	in item_list: 			line.append(clean_string(td["KE"]))
		if "PE" 	in item_list: 			line.append(clean_string(td["PE"]))
		if "Alpha" 	in item_list: 			line.append(clean_string(bd["Alpha"]))
		if "Beta" 	in item_list: 
			try: 
				line.append(clean_string(bd["Beta"]))
			except KeyError:
				line.append("-" + " "*14)
		if "N_MC" 	in item_list: 			line.append(r"$%15g$" % int(bd["N_MC"]))
		if "AcceptanceRate" in item_list: 	line.append(r"$%15.3g$" % float(td["AcceptanceRate"]))
		value_matrix.append(line)
	return value_matrix

def particle_info_getter(line):
	line = line.split("/")[-1]
	line = line.split("_")
	particle_number = [item for item in line if "Particle" in item]
	particle_omega = [item for item in line if "omega" in item]
	return [getDigit(particle_number[0]),getFraction(particle_omega[0])]

def sort_terminal_files(files):
	# Sorts after particle, then after omega
	new_list = []
	for file in files:
		new_list.append([file,particle_info_getter(file)])
	new_list.sort(key = lambda x: (x[1][0],-x[1][1]))
	for i in range(len(new_list)):
		files[i] = new_list[i][0]
	return files

# ===================================================================================================================================================
# Cases =============================================================================================================================================
# ===================================================================================================================================================

# Hard-coded ========================================================================================================================================
hc_folders = ["2e_JastrowWithCoulomb","2e_Jastrow","2e_plain"]
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
		line.append("-" + " "*14)
	line.append(r"$%15g$" % int(bd["N_MC"]))
	line.append(r"$%15.3g$" % float(td["AcceptanceRate"]))
	hc_value_matrix.append(line)

hc_tab = create_table(hc_value_matrix, row_values=hc_header_string, column_values=hc_col_string)
# print "Hard-coded table:\n"
# print hc_tab

# No interaction ====================================================================================================================================
no_int_folders = ["no_interaction"]
no_int_blocking_files = get_files(no_int_folders,blocking_output_directory)[0]
no_int_terminal_files = sort_terminal_files(get_files(no_int_folders,terminal_output_directory)[0])
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
no_int_tab = create_table(no_int_value_matrix, row_values=no_int_header_string)
# print "No interaction table"
# print no_int_tab

# No importance sampling ============================================================================================================================
no_imp_folders = ["no_imp"]
no_imp_blocking_files = get_files(no_imp_folders,blocking_output_directory)[0]
no_imp_terminal_files = sort_terminal_files(get_files(no_imp_folders,terminal_output_directory)[0])
no_imp_blocking_data = retrieve_file_contents(no_imp_blocking_files)
no_imp_terminal_data = retrieve_file_contents(no_imp_terminal_files,True)

# Parameter list
no_imp_param_list = []
for i in xrange(len(no_imp_blocking_data)/2):
	no_imp_param_list.append([	no_imp_blocking_data[2*i]["N"], 
							clean_omega(no_imp_blocking_data[2*i]["Omega"]), 
							clean_string(no_imp_blocking_data[2*i+1]["Alpha"]), # with Jastrow
							clean_string(no_imp_blocking_data[2*i+1]["Beta"]), # with Jastrow
							clean_string(no_imp_blocking_data[2*i]["Alpha"])]) # No Jastrow
no_imp_param_header = [r"$N$",
						r"$\omega$", 
						r"$\alpha$", 
						r"$\beta$", 
						r"$\alpha_\text{no Jastrow}$"]
no_imp_param_tab = create_table(no_imp_param_list,row_values=no_imp_param_header)
# print "Paramater table:"
# print no_imp_param_tab

# Energies list
new_no_imp_blocking_data = []
new_no_imp_terminal_data = []
for i in xrange(len(no_imp_blocking_data)): # Removing all data without jastrow factor
	if "Beta" in no_imp_blocking_data[i].keys():
		new_no_imp_blocking_data.append(no_imp_blocking_data[i])
		new_no_imp_blocking_data[-1]["EwoJ"] = no_imp_blocking_data[i-1]["Alpha"]
		new_no_imp_blocking_data[-1]["ERef"] = ERef[i/2]
		new_no_imp_terminal_data.append(no_imp_terminal_data[i])

no_imp_energy_header = [r"$N$", 
						r"$\omega$", 
						r"$\langle E \rangle$", 
						r"$\langle E \rangle_\text{no Jastrow}$", 
						r"$\langle E \rangle_1$", 
						r"$\sigma^2$", 
						r"$\langle K\rangle$", 
						r"$\langle V \rangle$", 
						r"$N_{MC}$",
						r"Acceptance"]

no_imp_value_matrix = create_value_matrix(new_no_imp_blocking_data,new_no_imp_terminal_data,["N","Omega","E","EwoJ","ERef","EVar","KE","PE","N_MC","AcceptanceRate"])
no_imp_tab = create_table(no_imp_value_matrix,row_values=no_imp_energy_header)
# print "No importance sampling table: "
# print no_imp_tab

# Importance sampling ===============================================================================================================================
#NOTE: Should create a matrix comparing this to regular sampling with energies with and without Jastrow, and variances especially
imp_folders = ["imp"]
imp_blocking_files = get_files(imp_folders,blocking_output_directory)[0]
imp_terminal_files = sort_terminal_files(get_files(imp_folders,terminal_output_directory)[0])
imp_blocking_data = retrieve_file_contents(imp_blocking_files)
imp_terminal_data = retrieve_file_contents(imp_terminal_files,True)
new_imp_data = []

imp_data = zip( imp_blocking_data, imp_terminal_data, no_imp_blocking_data[:10], no_imp_terminal_data[:10], ERef[:10])
imp_matrix = []
for i in xrange(len(imp_data)/2):
	new_dict = {}
	line = []
	new_dict["Omega"] = imp_data[2*i][0]["Omega"]
	line.append(clean_omega(imp_data[2*i][0]["Omega"]))
	new_dict["E"] = imp_data[2*i+1][0]["E"]
	line.append(clean_string(imp_data[2*i+1][0]["E"]))
	new_dict["EnoJ"] = imp_data[2*i][0]["E"]
	line.append(clean_string(imp_data[2*i][0]["E"]))
	new_dict["EnoImp"] = imp_data[2*i+1][2]["E"]
	line.append(clean_string(imp_data[2*i+1][2]["E"]))
	new_dict["EnoJnoImp"] = imp_data[2*i][2]["E"]
	line.append(clean_string(imp_data[2*i][2]["E"]))
	new_dict["ERef"] = ERef[i]
	line.append(clean_string(ERef[i]))
	new_dict["EVar"] = imp_data[2*i][0]["EVar"]
	line.append(clean_string(imp_data[2*i][0]["EVar"]))
	new_dict["EVarNoImp"] = imp_data[2*i][2]["EVar"]
	line.append(clean_string(imp_data[2*i][2]["EVar"]))
	new_dict["AcceptanceRate"] = imp_data[2*i][1]["AcceptanceRate"]
	line.append(clean_string("%15.3g" % float(imp_data[2*i][1]["AcceptanceRate"])))
	imp_matrix.append(line)
	new_imp_data.append(new_dict)

imp_header = [	r"$\omega$", 
				r"$\langle E \rangle$", 
				r"$\langle E \rangle_\text{no Jastrow}$", 
				r"$\langle E \rangle_\text{no imp.}$", 
				r"$\langle E \rangle_\text{no_Jastrow,no imp.}$", 
				r"$\langle E \rangle_1$", 
				r"$\sigma^2$", 
				r"$\sigma^2_\text{no imp.}$", 
				r"Imp. Acceptance"]

imp_tab = create_table(imp_matrix,row_values=imp_header)

###### FINAL RESULTS
print "Hard-coded table:\n"
print hc_tab
print "No interaction table"
print no_int_tab
print "Paramater table:"
print no_imp_param_tab
print "No importance sampling table: "
print no_imp_tab
print "Importance sampling versus regular matrix"
print imp_tab