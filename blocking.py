import numpy as np, matplotlib.pyplot as plt, os, sys, multiprocessing, re, time

def block(N_block_size):
	"""
	Function to be used in the parallel blocking method.
	Argument:
		N_block_sizes	: number of blocks we will split our data set into 
	"""
	print "Process %d now running for: %d" % (os.getpid(), N_block_size)
	block_size = N / N_block_size # Gets the size of each block
	temporary_average_values = np.zeros(N_block_size)
	for j in xrange(0,N_block_size):
		temporary_average_values[j] = np.sum(data[j*block_size:(j+1)*block_size])/block_size
	E = np.sum(temporary_average_values)/float(N_block_size)
	ESquared = np.sum(temporary_average_values**2)/float(N_block_size)
	del temporary_average_values # FOR RELEAVING MEMORY?
	# variance_value = (ESquared - E*E)/N_block_size
	return (ESquared - E*E)/N_block_size, block_size, N_block_size

num_processors = 4
verbose = True
dry_run = False # For just testing if everyting apart from blocking and plotting works
N_electron_values = [2, 6, 12, 20]
omega_values = [1.0, 0.5, 0.1, 0.05, 0.01]
res_b_sizes = [2e4, 3e5, 5e5, 5e5, 5e5]

def config_string(_folder,_n,_omega,_beta):
	return "folder: %s N_electrons: %2g Omega: %4g Beta: %2g" % (_folder, _n, _omega, _beta)

def configuration_error(error_message, _folder,_n,_omega,_beta):
	print "%s. Configuration does not exist: " % error_message + config_string(_folder,_n,_omega,_beta)

def factors(number):
    b = np.arange(1,number+1)
    res, = np.where((number % b) == 0)
    return np.array(res + 1)

for data_sub_folder in ["imp","no_imp","2e_plain","2e_jastrow","2e_jastrowWithCoulomb","no_interaction"]:
	folder_name = "output/" + data_sub_folder
	output_file = open("blocking_output/" + data_sub_folder + "/blocking_data_" + data_sub_folder + ".txt", "w")
	for N_electron in N_electron_values:
		for omega, resulting_block_size in zip(omega_values, res_b_sizes):
			for beta in [False,True]:

				pre_time = time.clock()
				figure_name = "electron%d_omega%f_beta%.4f" % (N_electron,omega,beta)
				try:
					raw_file_list = os.listdir(folder_name)
				except OSError:
					configuration_error("Warning: folder missing",data_sub_folder,N_electron,omega,beta)
					continue

				if len(raw_file_list) == 0:
					configuration_error("Warning: file list is empty", data_sub_folder,N_electron,omega,beta)
					continue

				# Building correct file list
				getDigit = lambda s : float(re.findall(r"[\d]+",s)[0])
				getFraction = lambda s : float(re.findall(r"[\d.]+",s)[0])
				file_list = []
				for file in raw_file_list:
					file_beta = None
					file_values = file.split('_')
					file_processor = getDigit(file_values[0][1:]) # the [1:] ensures we dont pick up any unnessecary numbers from folder
					file_particles = getDigit(file_values[1])
					file_MC = getDigit(file_values[2])
					file_omega = getFraction(file_values[3])
					file_alpha = getFraction(file_values[4])
					if beta:
						if "beta" in file:
							file_beta = getFraction(file_values[5])
						else:
							continue
					else:
						if "beta" in file:
							continue
					if file_particles != N_electron:
						continue
					if file_omega != omega:
						continue
					file_list.append(file)
				if len(file_list) == 0:
					if verbose: configuration_error("No data files found", data_sub_folder,N_electron,omega,beta)
					continue

				data_files = []
				# Constants to be written to file++
				alpha_value = 0
				beta_value = 0
				N_MC = 0

				# Retrieves data
				if verbose: print "\nLoading data for configuration: %s" % config_string(folder_name,N_electron,omega,beta)
				for file in file_list:
					file_beta = None
					file_values = file.split('_')
					file_processor = getDigit(file_values[0][1:]) # the [1:] ensures we dont pick up any unnessecary numbers from folder
					file_particles = getDigit(file_values[1])
					file_MC = getDigit(file_values[2])
					file_omega = getFraction(file_values[3])
					file_alpha = getFraction(file_values[4])
					if beta:
						file_beta = getFraction(file_values[5])
					data_files.append(np.fromfile(folder_name + "/" + file))
					alpha_value = file_alpha
					N_MC = file_MC
					print_string = "Data loaded from processor %d: N = %2d MC Cycles = %10d Omega = %10.8f Alpha = %10.8f" % (file_processor, file_particles, file_MC, file_omega, file_alpha)
					if beta and (file_beta != None):
						beta_value = file_beta
						print_string += " Beta = %10.8f" % file_beta 
					if verbose: print print_string

				if len(data_files) == 0:
					if verbose: configuration_error("No data files found. SHOULD NOT SEE THIS MESSAGE", data_sub_folder,N_electron,omega,beta)
					continue

				data = np.concatenate(data_files)
				del data_files
				print_string = "All data loaded from %s for %3d electrons, N_MC = %5g, omega = %8.6f, alpha = %10.8f" % (folder_name, N_electron, N_MC, omega, alpha_value)
				if beta: print_string += " Beta = %10.8f" % beta_value
				if verbose: print print_string

				# Setting up blocks ---------------------------------------
				N = len(data)
				N_block_sizes = [] # Is the number of blocks we divide into
				temp_N = N
				while temp_N >= 1e6:
					temp_N /= 10
				# for i in xrange(1,temp_N/2): # Finds the number of blocks to divide the data into.
				# 	if (N % i == 0):
				# 		N_block_sizes.append(N/i)
				N_block_sizes = factors(temp_N)[::-1]
				N_block_sizes = N_block_sizes[:-2] # Skipping first four and last one(since the last one varies alot and is unreliable)
				N_blocks = len(N_block_sizes)
				print N_block_sizes
				if dry_run: continue

				# Run parallelized blocking -------------------------------
				if N >= 1e8: num_processors = 2 # In order to not kill my computer
				if verbose: "Performing blocking with %d processors" % num_processors
				pool = multiprocessing.Pool(processes=num_processors)
				res = pool.map(block,N_block_sizes)
				pool.close()
				res = np.asarray(res)
				N_block_sizes = res[:,2]
				block_sizes = res[:,1]
				variance_values = res[:,0]

				# Finding the energy mean ---------------------------------
				energy_mean = np.mean(data)
				energy_variance = block(int(N/resulting_block_size))[0]
				data_string = "N %4d N_MC %10d E %18.16f EVar %12.g EStd %12.g Omega %6.4f Alpha %12.10f" % (N_electron, N_MC, energy_mean, energy_variance, np.sqrt(energy_variance), omega, alpha_value)
				if beta: data_string += " Beta %.10f" % beta_value
				output_file.write(data_string + "\n")

				# Plotting ------------------------------------------------
				plt.clf()
				plt.cla()
				plt.semilogx(block_sizes[::-1], variance_values[::-1])
				plt.xlabel("Block size")
				plt.ylabel(r"Variance $\sigma^2$")
				plt.title(r"Blocking for $%d$ electrons, $\omega=%.2f$, $N_{MC}=%d$" % (N_electron, omega, N_MC))
				plt.grid(True)
				plt.savefig("figures/%s/%s.png" % (data_sub_folder, figure_name),dpi=300)
				plt.close()
				if verbose: print "figures/%s/%s.png plotted" % (data_sub_folder, figure_name)
				# plt.show()

				del res, N_block_sizes, block_sizes, variance_values, data
				post_time = time.clock()
				print "Blocking complete. Time taken: %s seconds" % (post_time - pre_time)
	output_file.close()

# # Settings ------------------------------
# num_processors = 4
# n_electrons = 2
# figure_name = "2electron"
# folder_name = "output_2e_MC1e7/"

# file_list = os.listdir("output_2el_run2/") # For Ubuntu
# file_list = os.listdir(folder_name) # For Mac
# # Loading files -------------------------------------------
# data_files = []
# for file in file_list:
# 	data_files.append(np.fromfile(folder_name + file))
# data = np.concatenate(data_files)
# print "Data loaded from %s" % folder_name

# # Setting up blocks ---------------------------------------
# N = len(data)
# N_block_sizes = [] # Is the number of blocks we divide into
# for i in xrange(1,N/2): # Finds the number of blocks to divide the data into.
# 	if (N % i == 0):
# 		N_block_sizes.append(N/i)
# # REMOVE START BLOCKS AS THEY TAKE LONG TIME?
# N_blocks = len(N_block_sizes)

# # Non-parallelized version ------------------------------
# block_size_values = np.zeros(N_blocks)
# variance_values = np.zeros(N_blocks)
# for i in xrange(N_blocks):
# 	block_size = N / N_block_sizes[i]
# 	print i, block_size, N_block_sizes[i]
# 	temporary_average_values = np.zeros(N_block_sizes[i])
# 	for j in xrange(0,N_block_sizes[i]):
# 		temporary_average_values[j] = np.sum(data[j*block_size:(j+1)*block_size])/block_size
# 	E = np.sum(temporary_average_values)/float(N_block_sizes[i])
# 	ESquared = np.sum(temporary_average_values**2)/float(N_block_sizes[i])
# 	variance_values[i] = (ESquared - E*E)/N_block_sizes[i]
# 	block_size_values[i] = block_size

# # Run parallelized blocking -------------------------------
# pool = multiprocessing.Pool(processes=num_processors)
# res = pool.map(block,N_block_sizes)
# pool.close()
# res = np.asarray(res)
# N_block_sizes = res[:,2]
# block_sizes = res[:,1]
# variance_values = res[:,0]

# # Finding the energy mean ---------------------------------
# energy_mean = np.mean(data)
# print "E = ", energy_mean

# # Plotting ------------------------------------------------
# plt.semilogx(block_sizes[::-1], variance_values[::-1])
# plt.xlabel("Block size")
# plt.ylabel(r"Variance $\sigma$")
# plt.title(r"Blocking size versus variance for %d electrons" % n_electrons)
# plt.grid(True)
# plt.savefig("figures/%s.png" % figure_name,dpi=300)
# # plt.show()