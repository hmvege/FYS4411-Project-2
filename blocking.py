import numpy as np, matplotlib.pyplot as plt, os, sys

file_list = os.listdir("output_2el_run2/")
data_files = []
for file in file_list:
	data_files.append(np.fromfile("output_2el_run2/" + file))

data = np.concatenate(data_files)

N = len(data)
N_blocks = 200
min_block_size = 100
max_block_size = len(data)/100
block_step = int((max_block_size - min_block_size)/N_blocks + 1)

print np.var(data)/float(N)
print np.mean(data)
print """
Data set size: 		%d
Min block size: 	%d
Max block size: 	%d
Block step:		%d
Number of blocks: 	%d""" % (N, min_block_size, max_block_size, 
	N_blocks, block_step)

def block_mean(x, n_start, n_stop):
	return np.mean(x[n_start:n_stop])

def block(x, n, b_size):
	num_blocks = int(n/b_size)
	block_values = np.zeros(num_blocks)
	for i in xrange(num_blocks-1):
		block_values[i] = np.mean(x[i*b_size:(i+1)*b_size])
	block_average = np.mean(block_values)
	# block_variance = np.var(block_values)
	block_variance = np.sum(block_values*block_values)/len(block_values) - block_average*block_average
	return block_average, block_variance

average_values = np.zeros(N_blocks)
variance_values = np.zeros(N_blocks)
block_sizes = np.zeros(N_blocks)

for i in xrange(N_blocks):
	block_size = min_block_size + block_step*i
	average_values[i], variance_values[i] = block(data, N, block_size)
	variance_values[i] = np.sqrt(variance_values[i]/(N/float(block_size) - 1.0))
	block_sizes[i] = block_size
print "E = ", np.mean(average_values)
plt.plot(block_sizes, variance_values)
plt.show()