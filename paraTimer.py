import numpy as np

para = np.array([22.313027,22.484323,22.220841,22.757712])
non_para = np.array([75.566675,74.749889,75.247932,75.280813])
processors = 4

def speedup(T1, Tp):
	mean_T1 = np.mean(T1)
	mean_Tp = np.mean(Tp)
	delta_T1 = np.std(T1)
	delta_Tp = np.std(Tp)
	delta_S = np.sqrt( (delta_T1/mean_T1)**2 + (delta_Tp/mean_Tp)**2 )
	S = mean_T1/mean_Tp
	return S, delta_S

def f(T1, Tp, p):
	S = np.mean(T1)/np.mean(Tp)
	return (1-S)/(S*(1/p-1))

def delta_f(T1, Tp, p):
	# mean_T1 = np.mean(T1)
	# mean_Tp = np.mean(Tp)
	# S = mean_T1/mean_Tp
	# delta_T1 = np.std(T1)
	# delta_Tp = np.std(Tp)
	# delta_S = np.sqrt( (delta_T1/mean_T1)**2 + (delta_Tp/mean_Tp)**2 )
	S, delta_S = speedup(T1, Tp)
	return f(T1, Tp, p)*delta_S / (S*(S-1))

speed_up, speed_up_delta = speedup(non_para,para)

print "Speed-up: T1/Tp = %g +/- %g" % (speed_up,speed_up_delta)
print "Fraction that is parallelizable: f = %g +/- %g" % (
	f(non_para,para,processors),delta_f(non_para,para,processors))

F = (1 - (non_para/para))/((non_para/para)*(1/processors-1))
print np.mean(non_para/para), np.std(non_para/para)
print np.mean(F), np.std(F)