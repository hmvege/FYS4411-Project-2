import numpy as np, matplotlib.pyplot as plt, os, sys, multiprocessing

# Settings ------------------------------
num_processors = 4
n_electrons = 2
figure_name = "2electron"
folder_name = "output_2e_MC1e7/"
# file_list = os.listdir("output_2el_run2/") # For Ubuntu
file_list = os.listdir(folder_name) # For Mac

# Loading files -------------------------------------------
data_files = []
for file in file_list:
	data_files.append(np.fromfile(folder_name + file))
data = np.concatenate(data_files)
print "Data loaded from %s" % folder_name

# Setting up blocks ---------------------------------------
N = len(data)
N_block_sizes = [] # Is the number of blocks we divide into
for i in xrange(1,N/2): # Finds the number of blocks to divide the data into.
	if (N % i == 0):
		N_block_sizes.append(N/i)
N_blocks = len(N_block_sizes)

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

def block(N_block_size):
	"""
	Function to be used in the parallel blocking method.
	Argument:
		N_block_sizes	: number of blocks we will split our data set into 
	"""
	print "Now running for: ", N_block_size
	block_size = N / N_block_size # Gets the size of each block
	temporary_average_values = np.zeros(N_block_size)
	for j in xrange(0,N_block_size):
		temporary_average_values[j] = np.sum(data[j*block_size:(j+1)*block_size])/block_size
	E = np.sum(temporary_average_values)/float(N_block_size)
	ESquared = np.sum(temporary_average_values**2)/float(N_block_size)
	del temporary_average_values # FOR RELEAVING MEMORY?
	# variance_value = (ESquared - E*E)/N_block_size
	return (ESquared - E*E)/N_block_size, block_size, N_block_size

# Run parallelized blocking -------------------------------
pool = multiprocessing.Pool(processes=num_processors)
res = pool.map(block,N_block_sizes)
res = np.asarray(res)
N_block_sizes = res[:,2]
block_sizes = res[:,1]
variance_values = res[:,0]

# Finding the energy mean ---------------------------------
energy_mean = np.mean(data)
print "E = ", energy_mean

# Plotting ------------------------------------------------
plt.semilogx(block_sizes[::-1], variance_values[::-1])
plt.xlabel("Block size")
plt.ylabel(r"Variance $\sigma$")
plt.title(r"Blocking size versus variance for %d electrons" % n_electrons)
plt.grid(True)
plt.savefig("figures/%s.png" % figure_name,dpi=300)
# plt.show()