import numpy as np, matplotlib.pyplot as plt, os, sys, multiprocessing, re, time

def block(N_block_size):
	"""
	Function to be used in the parallel blocking method.
	Argument:
		N_block_sizes	: number of blocks we will split our data set into 
	"""
	# print "Process %d now running for: %d" % (os.getpid(), N_block_size)
	block_size = N / N_block_size # Gets the size of each block
	temporary_average_values = np.zeros(N_block_size)
	for j in xrange(0,N_block_size):
		temporary_average_values[j] = np.sum(data[j*block_size:(j+1)*block_size])/block_size
	E = np.sum(temporary_average_values)/float(N_block_size)
	ESquared = np.sum(temporary_average_values**2)/float(N_block_size)
	del temporary_average_values # FOR RELEAVING MEMORY?
	# variance_value = (ESquared - E*E)/N_block_size
	return (ESquared - E*E)/N_block_size, block_size, N_block_size

N_electron_values = [2, 6, 12, 20]
omega_values = [1.0, 0.5, 0.1, 0.05, 0.01]
folder_name = "output_no_imp_run"
beta_vals = [True, False]
num_processors = 8

for N_electron in N_electron_values:
	for omega in omega_values:
		for beta in beta_vals:
			pre_time = time.clock()

			figure_name = "electron%d_omega%f_beta%.4f" % (N_electron,omega,beta)
			file_list = os.listdir(folder_name)
			data_files = []
			getDigit = lambda s : float(re.findall(r"[\d]+",s)[0])
			getFraction = lambda s : float(re.findall(r"[\d.]+",s)[0])
			for file in file_list:
				file_beta = None
				file_values = file.split('_')
				file_processor = getDigit(file_values[0])
				file_particles = getDigit(file_values[1])
				file_MC = getDigit(file_values[2])
				file_omega = getFraction(file_values[3])
				file_alpha = getFraction(file_values[4])
				if len(file_values) > 5:
					file_beta = getFraction(file_values[5])
				if ((N_electron == file_particles) and (file_omega == omega) and (((not beta) and (file_beta == None)) or (beta and (file_beta != None)))):
					data_files.append(np.fromfile(folder_name + "/" + file))
					N_MC = file_MC
					print_string = "Data loaded from processor %2d: N = %2d MC Cycles = %10d Omega = %10.8f Alpha = %10.8f" % (file_processor, file_particles, file_MC, file_omega, file_alpha)
					if (file_beta != None):
						print_string += " Beta = %10.8f" % file_beta 
					print print_string
			data = np.concatenate(data_files)
			print "Data loaded from %s" % folder_name

			# Setting up blocks ---------------------------------------
			N = len(data)
			N_block_sizes = [] # Is the number of blocks we divide into
			for i in xrange(1,N/2): # Finds the number of blocks to divide the data into.
				if (N % i == 0):
					N_block_sizes.append(N/i)
			N_block_sizes = N_block_sizes[8:-1] # Skipping first four and last one
			N_blocks = len(N_block_sizes)

			# Run parallelized blocking -------------------------------
			pool = multiprocessing.Pool(processes=num_processors)
			res = pool.map(block,N_block_sizes)
			pool.close()
			res = np.asarray(res)
			N_block_sizes = res[:,2]
			block_sizes = res[:,1]
			variance_values = res[:,0]

			# Finding the energy mean ---------------------------------
			energy_mean = np.mean(data)
			energy_variance = block(int(N/1e4))[0]
			print "E = %.10f +/- %g" % (energy_mean, np.sqrt(energy_variance))

			# Plotting ------------------------------------------------
			plt.clf()
			plt.cla()
			plt.semilogx(block_sizes[::-1], variance_values[::-1])
			plt.xlabel("Block size")
			plt.ylabel(r"Variance $\sigma^2$")
			plt.title(r"Blocking for $%d$ electrons, $\omega=%.2f$, $N_{MC}=%d$" % (N_electron, omega, N_MC))
			plt.grid(True)
			plt.savefig("figures/%s.png" % figure_name,dpi=300)
			plt.close()
			# plt.show()

			del res, N_block_sizes, block_sizes, variance_values, data
			post_time = time.clock()
			print "Blocking complete. Time taken: %s seconds" % (post_time - pre_time)


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