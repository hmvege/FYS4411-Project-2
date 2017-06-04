import numpy as np, matplotlib.pyplot as plt, os

file_list = os.listdir("output_2el_run2/")
data_files = []
for file in file_list:
	data_files.append(np.fromfile("output_2el_run2/" + file))

data = np.concatenate(data_files)

N = len(data)
print N
N_blocks = 200
min_block_size = 100
max_block_size = len(data)/100
block_step = int((max_block_size-min_block_size + 1)/N_blocks)

def block_mean(x, n_start, n_stop):
	return np.mean(x[n_start:n_stop])

def block(x, n, b_size):
	num_blocks = int(n/b_size)
	block_values = np.zeros(num_blocks)
	for i in xrange(num_blocks-1):
		block_values[i] = np.mean(x[i*b_size:(i+1)*b_size])
	block_average = np.mean(block_values)
	block_variance = np.var(block_values)/float(len(block_values))
	return block_average, block_variance

average_values = np.zeros(N_blocks)
variance_values = np.zeros(N_blocks)
block_sizes = np.zeros(N_blocks)

for i in xrange(N_blocks):
	block_size = min_block_size + block_step*i
	average_values[i], variance_values[i] = block(data, N, block_size)
	block_sizes[i] = block_size

plt.plot(block_sizes, variance_values)
plt.show()