#rewrite data import as 3d array
#rewrite derivatives using np.diff


import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt
import subprocess
import os
import matplotlib as mpl

font = {'family':'serif','weight':'normal','size':12}
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['axes.labelsize'] = 16
mpl.rc('font', **font)

#Standard values for [wb,wc,wa,h,ns]
central_values = [0.02238,0.125,0.0000001,67.7,0.96]
#Fine tuned step sizes to ensure well-behaved derivatives
fin_epsilon_definite = [0.00125,0.0051,0.001,0.08,0.025]
#Converts the step sizes to fractional step sizes
fin_epsilon = [fin_epsilon_definite[i]/central_values[i] for i in range(len(fin_epsilon_definite))]
#Creates upper and lower values for each parameter to obtain accurate derivatives
fin_varied = {key: [] for key in range(len(central_values))}
for key in fin_varied:
	#wa uses a one-sided derivative to avoid negative values
	if key == 2:
		fin_varied[key].append([central_values[key],central_values[key]*(1.0+fin_epsilon[key])])
	#All other parameters use two-sided derivatives
	else:
		fin_varied[key].append([central_values[key]*(1.0-fin_epsilon[key]),central_values[key]*(1.0+fin_epsilon[key])])
#Creates a list of masses even space in log-space
masses = np.logspace(-24,-30,100)

new_fin_epsilon = [0.05585344057193923,0.026041666666666668]
new_central_values = [0.02238,0.96]

a = np.arange(0,100,1)
b = np.arange(1,201,2)
c = np.arange(2,201,2)

epsilon_matrix = np.asarray(np.repeat(new_fin_epsilon,100)).reshape(2,100,1)
param_matrix = np.asarray(np.repeat(new_central_values,100)).reshape(2,100,1)
print epsilon_matrix
print np.shape(epsilon_matrix)

k_low = []
Pk_low = []
for i in range(len(b)):
	k_low.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/wb%s_m%s.dat' %(b[i],a[i]), usecols = [0]))
	Pk_low.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/wb%s_m%s.dat' %(b[i],a[i]), usecols = [1]))
for i in range(len(b)):
	Pk_low.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/ns%s_m%s.dat' %(b[i],a[i]), usecols = [1]))
k_low = np.asarray(k_low)
Pk_low =np.asarray(Pk_low).reshape(2,100,551)

k_high = []
Pk_high = []
for i in range(len(c)):
        k_high.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/wb%s_m%s.dat' %(c[i],a[i]), usecols = [0]))
        Pk_high.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/wb%s_m%s.dat' %(c[i],a[i]), usecols = [1]))
for i in range(len(c)):
	k_high.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/ns%s_m%s.dat' %(c[i],a[i]), usecols = [0]))
        Pk_high.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/ns%s_m%s.dat' %(c[i],a[i]), usecols = [1]))
k_high =np.asarray(k_high).reshape(2,100,551)
Pk_high =np.asarray(Pk_high).reshape(2,100,551)

def fin_deriv(Pk_high,Pk_low,central_value,step_size):
	param_deriv = (Pk_high-Pk_low)/(central_value*step_size)
	return param_deriv

deriv = fin_deriv(Pk_high,Pk_low,param_matrix,epsilon_matrix)
print np.shape(deriv)
plt.plot(k_high[1][0],deriv[1][0])
plt.xlim(0.0001,1.0)
plt.xscale('log')
plt.savefig('test_deriv.png')
subprocess.call('open test_deriv.png',shell = True)
"""
even = np.arange(0,len(Pk),2)
odd = np.arange(1,len(Pk)+1,2)

#test = [[np.subtract(Pk[i][even[j]],Pk[i][odd[j]]) for j in range(len(Pk[i]))] for i in range(len(Pk))]
#for i in range(len(Pk)):
#print test

#Initializes finite difference power spectrum dictionaries
k_wb = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
Pk_wb = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
k_wc = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
Pk_wc = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
k_wa = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
Pk_wa = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
k_ns = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
Pk_ns = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
k_h = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
Pk_h = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
k = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}
Pk = {key:{i:[] for i in range(2*len(fin_varied[1]))} for key in range(len(masses))}

#Fills the finite difference power spectrum dictionaries
for key in k:
        counter = 0
        for i in range(2*key+1,2*key+3):
                k_wb[key][counter],Pk_wb[key][counter] = loadtxt('mp/wb%s_m%s.dat' %(i,key), unpack = True, usecols=[0,1])
                k_wc[key][counter],Pk_wc[key][counter] = loadtxt('mp/wc%s_m%s.dat' %(i,key), unpack = True, usecols=[0,1])
                k_wa[key][counter],Pk_wa[key][counter] = loadtxt('mp/wa%s_m%s.dat' %(i,key), unpack = True, usecols=[0,1])
                k_ns[key][counter],Pk_ns[key][counter] = loadtxt('mp/ns%s_m%s.dat' %(i,key), unpack = True, usecols=[0,1])
                k_h[key][counter],Pk_h[key][counter] = loadtxt('mp/h%s_m%s.dat' %(i,key), unpack = True, usecols=[0,1])
                k[key][counter],Pk[key][counter] = loadtxt('mp/m%s.dat' %(key), unpack = True, usecols=[0,1])
                counter += 1


#Survey parameters
n = 0.0001
kmax = 0.11
#kmax = 0.2
V_survey = 1.0
kmin = 2*math.pi/(((V_survey)**(1/3))*1000)
#Limits the wavenumbers to kmin < k < kmax
k_in_range = {key:[k[key][1][i] for i in range(len(k[key][1])) if kmin < k[key][1][i] < kmax] for key in k}
#Creates evenly spaced wavenumbers in the range kmin < k < kmax
even_k_grid = {key:np.linspace(kmin,kmax,len(k_in_range[key])) for key in k}
#Selects the Pk's that correspond to the relevant k's
Pk_in_range = {key:[Pk[key][1][i] for i in range(len(k[key][1])) if kmin < k[key][1][i] < kmax] for key in k_wb}
#Interpolates the Pk's to occur at the evenly spaced k's specified in even_k_grid
Pk_in_range = {key:np.interp(even_k_grid[key],k_in_range[key],Pk_in_range[key]) for key in k_wb}
#Computes the effective survey volume given by eq.10 of Seo & Eisenstein 2003
V_eff = {key:[np.square(((n*Pk_in_range[key][i])/(n*Pk_in_range[key][i]+1)))*((1000**3)*V_survey) for i in range(len(Pk_in_range[key]))] for key in k_wb}
#Selects the parameter derivative values that occure within kmin < k < kmax, then interpolates to match with  even_k_grid k's
wb_deriv_in_range = {key:[wb_deriv[key][i] for i in range(len(wb_deriv[key])) if kmin < k_wb[key][1][i] < kmax] for key in k_wb}
wb_deriv_in_range = {key:np.interp(even_k_grid[key],k_in_range[key],wb_deriv_in_range[key]) for key in k_wb}
wc_deriv_in_range = {key:[wc_deriv[key][i] for i in range(len(wc_deriv[key])) if kmin < k_wc[key][1][i] < kmax] for key in k_wb}
wc_deriv_in_range = {key:np.interp(even_k_grid[key],k_in_range[key],wc_deriv_in_range[key]) for key in k_wb}
wa_deriv_in_range = {key:[wa_deriv[key][i] for i in range(len(wa_deriv[key])) if kmin < k_wa[key][1][i] < kmax] for key in k_wb}
wa_deriv_in_range = {key:np.interp(even_k_grid[key],k_in_range[key],wa_deriv_in_range[key]) for key in k_wb}
ns_deriv_in_range = {key:[ns_deriv[key][i] for i in range(len(ns_deriv[key])) if kmin < k_ns[key][1][i] < kmax] for key in k_wb}
ns_deriv_in_range = {key:np.interp(even_k_grid[key],k_in_range[key],ns_deriv_in_range[key]) for key in k_wb}
h_deriv_in_range = {key:[h_deriv[key][i] for i in range(len(h_deriv[key])) if kmin < k_h[key][1][i] < kmax] for key in k_wb}
h_deriv_in_range = {key:np.interp(even_k_grid[key],k_in_range[key],h_deriv_in_range[key]) for key in k_wb}
#Computes the uniform step size between k's
dk = {key:(kmax-kmin)/(len(k_in_range[key])-1) for key in k_wb}
#Initializes an empty matrix for each mass that will become the fisher matrices
matrix_dict = {key:np.zeros((5,5)) for key in k_wb}
#Initializes a dictionary whose keys are mass runs, then each mass has empty lists to be filled by each parameter
in_range_dict = {key:{i:[] for i in range(len(central_values))} for key in k_wb}

#Fills in_range_dict -- key loops through the masses, and then each mass has 5 sub-entries, one for each parameter
for key in in_range_dict:
	in_range_dict[key][0] = wb_deriv_in_range[key]
	in_range_dict[key][1] = wc_deriv_in_range[key]
	in_range_dict[key][2] = wa_deriv_in_range[key]
	in_range_dict[key][3] = ns_deriv_in_range[key]
	in_range_dict[key][4] = h_deriv_in_range[key]

#Uses the entries in in_range_dict to fill the fisher matrix
for key in in_range_dict:
	for m in range(len(in_range_dict[1])):
		for n in range(len(in_range_dict[1])):
			matrix_dict[key][m][n] = np.trapz([in_range_dict[key][m][i]*in_range_dict[key][n][i]*V_eff[key][i]*np.square(k_in_range[key][i])/(4*np.square(math.pi)*np.square(Pk_in_range[key][i])) for i in range(len(k_in_range[key]))],k_in_range[key],dk[key])

#Cholesky decomposition to invert the fisher matrices (matrix_dict) and obtain the error matrices
#Cannot perform a simple np.invert because the fisher matrices are ill-conditioned
inverted_dict = {key:np.dot(np.linalg.inv(np.linalg.cholesky(matrix_dict[key]).transpose()),np.linalg.inv(np.linalg.cholesky(matrix_dict[key]))) for key in k_wb}

#Dictionary of the variance of each parameter: key for each parameter, then a mass run entries for each
variance_dict = {i+1:[np.sqrt(inverted_dict[key][i][i]) for key in k_wb] for i in range(len(central_values))}
wb_variance = [np.sqrt(inverted_dict[key][0][0]) for key in k_wb]
wc_variance = [np.sqrt(inverted_dict[key][1][1]) for key in k_wb]
wa_variance = [np.sqrt(inverted_dict[key][2][2]) for key in k_wb]
ns_variance = [np.sqrt(inverted_dict[key][3][3]) for key in k_wb]
h_variance = [np.sqrt(inverted_dict[key][4][4]) for key in k_wb]

fig,ax = plt.subplots()
ax.plot(k_wb[50][0],wb_deriv[50],label = 'wb')
ax.plot(k_wc[50][0],wc_deriv[50],label = 'wc')
#ax.plot(k_wa[50][0],wa_deriv[50],label = 'wa')
ax.plot(k_ns[50][0],ns_deriv[50],label = 'ns')
ax.plot(k_h[50][0],h_deriv[50],label = 'h')
ax.set_xlim(0.0001,1.0)
ax.set_xscale('log')
plt.legend(frameon = False, loc = 'lower left')
plt.savefig('test.png')
subprocess.call('open test.png',shell = True)
plt.close('all')
"""
