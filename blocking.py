import numpy as np, matplotlib.pyplot as plt, os

file_list = os.listdir("output/")
print file_list


data = np.fromfile("output/NElectron_Particle6_MC100000_omega1.000000_alpha1.000000_beta0.400000")

print len(data)
print np.mean(data)