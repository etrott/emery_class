#initial fisher code  swaps the ns and h denominators in the derivative step (fixed here)

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

def load(param_string,param_index_list,mass_index_list):
	k = []
	Pk = []
	for i in range(len(param_index_list)):
		k.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/%s%s_m%s.dat' %(param_string,param_index_list[i],mass_index_list[i]), usecols = [0]))
		Pk.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/%s%s_m%s.dat' %(param_string,param_index_list[i],mass_index_list[i]), usecols = [1]))
	return k,Pk

k_low = []
Pk_low = []
for i in range(len(b)):
	k_low.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/wb%s_m%s.dat' %(b[i],a[i]), usecols = [0]))
	Pk_low.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/wb%s_m%s.dat' %(b[i],a[i]), usecols = [1]))
"""
k_low,Pk_low = load('wb',b,a)
"""
for i in range(len(b)):
        k_low.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/ns%s_m%s.dat' %(b[i],a[i]), usecols = [0]))
	Pk_low.append(loadtxt('/Users/etrott12/Dropbox/emery/axionCAMB/mp/ns%s_m%s.dat' %(b[i],a[i]), usecols = [1]))
k_low = np.asarray(k_low).reshape(2,100,551)
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

#Survey parameters
n = np.asarray(0.0001).reshape(1,1,1)
kmax = 0.11
V_survey = np.asarray(1.0).reshape(1,1,1)
kmin = 2*np.pi/(((V_survey)**(1/3))*1000)
#Limits the wavenumbers to kmin < k < kmax
bool = np.logical_and(np.greater(k_low,kmin),np.less(k_low,kmax))
k_ir = k_low[bool].reshape(2,100,143)
Pk_ir = Pk_low[bool].reshape(2,100,143)
deriv_ir = deriv[bool].reshape(2,100,143)

"""
k_ir = k_low[np.all([k_low > kmin, k_low < kmax])]
deriv_ir = deriv[np.all([k_low > kmin, k_low < kmax])]
even_k = np.linspace(k_ir[0][0][0],k_ir[0][0][-1],len(k_low[0][0]))
deriv_int = np.interp(even_k,k_ir,deriv_ir)
test = []
for i in range(len(deriv)):
	for j in range(len(deriv[i])):
		test.append(np.interp(even_k,k_ir[i][j],deriv_ir[i][j]))
test = np.asarray(test).reshape(2,100,551)

plt.plot(k_ir[1][0],deriv_ir[1][0])
plt.xlim(0.0001,1.0)
plt.xscale('log')
plt.savefig('test_deriv.png')
subprocess.call('open test_deriv.png',shell = True)
"""

#Computes the effective survey volume given by eq.10 of Seo & Eisenstein 2003
V_eff = np.square((n*Pk_ir)/(n*Pk_ir+1))*V_survey
dk = k_ir[0][0][1]-k_ir[0][0][0]

k_matrix = k_ir[0].reshape(1,100,143)

matrix = []
for i in range(len(deriv_ir)):
	for j in range(len(deriv_ir)):
		matrix.append(np.trapz(deriv_ir[i]*deriv_ir[j]*V_eff[i]*np.square(k_ir[0])/(4*np.pi**2),k_matrix,dk))
matrix = np.transpose(matrix).reshape(100,2,2)
invert = []
for i in range(100):
	invert.append(np.dot(np.linalg.inv(np.linalg.cholesky(matrix[i]).transpose()),np.linalg.inv(np.linalg.cholesky(matrix[i]))))
new = np.reshape(invert,(100,4)).transpose()

plt.plot(masses,np.sqrt(new[0]),color = 'r',marker = '_',label = r'$\omega_b$')
plt.plot(masses,np.sqrt(new[3]),color ='b',marker = '_',label= r'$n_s$')
plt.fill_between(masses,np.sqrt(new[0]),alpha = 0.4,color = 'r')
plt.fill_between(masses,np.sqrt(new[3]),alpha = 0.4,color = 'b')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'mass (eV)')
plt.ylabel(r'$1\sigma$')
plt.title('Forecasted Constraints')
plt.legend(frameon = False, loc = 'upper right')
plt.savefig('test_deriv.png')
subprocess.call('open test_deriv.png',shell = True)
