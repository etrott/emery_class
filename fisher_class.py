
import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt
import subprocess
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

def run_CLASS(param,value):
    cosmo = Class()
    params = {'output':'tCl,pCl,lCl','lensing':'yes',param:value}
    cosmo.set(params)
    cosmo.compute()
    output = cosmo.lensed_cl(2500)
    return output['tt']

param_list = ['Omega_b','Omega_cdm','H0']
central_values = [0.022032,0.12038,67.556]
#Fine tuned step sizes to ensure well-behaved derivatives
fin_epsilon_definite = [0.00125,0.005,0.08]
#Converts the step sizes to fractional step sizes
fin_epsilon = [fin_epsilon_definite[i]/central_values[i] for i in range(len(fin_epsilon_definite))]
#Creates a list of masses even space in log-space
#masses = np.logspace(-24,-30,10)
masses = [1.0]
ell = np.arange(0,2501).reshape(1,2501)
param_matrix = np.asarray(central_values).reshape(len(param_list),1)
epsilon_matrix = np.asarray(fin_epsilon).reshape(len(param_list),1)
param_values = []
for i in range(len(param_list)):
    param_values.append((central_values[i]*(1.0-fin_epsilon[i]),central_values[i]*(1.0+fin_epsilon[i])))
param_values = np.reshape(param_values,(len(param_list),2))

"""
#runs CLASS and saves the outputs to npy files
cl_low = []
cl_high = []
for i in range(len(param_list)):
    cl_low.append(run_CLASS(param_list[i],param_values[i][0]))
    cl_high.append(run_CLASS(param_list[i],param_values[i][1]))
cl_low = np.reshape(cl_low,(len(param_list),2501))
cl_high = np.reshape(cl_high,(len(param_list),2501))
np.save('cl_low',cl_low)
np.save('cl_high',cl_high)
"""

#imports CLASS data
cl_low = np.load('cl_low.npy')
cl_high = np.load('cl_high.npy')

cl_low = np.reshape(cl_low,(len(fin_epsilon),len(masses),len(ell[0])))
cl_high = np.reshape(cl_high,(len(fin_epsilon),len(masses),len(ell[0])))

def fin_deriv(cl_high,cl_low,central_value,step_size):
	param_deriv = (cl_high-cl_low)/(central_value*step_size)
	return param_deriv
deriv = fin_deriv(cl_high,cl_low,param_matrix,epsilon_matrix)

n = np.asarray(0.0001).reshape(1,1,1)
kmax = 0.11
V_survey = np.asarray(1.0).reshape(1,1,1)
kmin = 2*np.pi/(((V_survey)**(1/3))*1000)
"""
bool = np.logical_and(np.greater(k_low,kmin),np.less(k_low,kmax))
k_ir = k_low[bool].reshape(len(fin_epsilon),len(masses),143)
Pk_ir = Pk_low[bool].reshape(len(fin_epsilon),len(masses),143)
deriv_ir = deriv[bool].reshape(len(fin_epsilon),len(masses),143)
"""


#Computes the effective survey volume given by eq.10 of Seo & Eisenstein 2003
V_eff = np.square((n*cl_low)/(n*cl_low+1))*V_survey
dk = ell[0][1]-ell[0][0]
ell_matrix = ell[0].reshape(1,len(masses),len(ell[0]))
matrix = []
for i in range(len(deriv)):
	for j in range(len(deriv)):
		matrix.append(np.trapz(deriv[i]*deriv[j]*V_eff[i]*np.square(ell[0])/(4*np.pi**2),ell_matrix,dk))
matrix = np.transpose(matrix)
matrix = np.reshape(matrix[0][0],(len(param_list),len(param_list)))

#matrix = np.transpose(matrix).reshape(len(masses),len(fin_epsilon),len(fin_epsilon))

#invert = []
#for i in range(100):
#	invert.append(np.dot(np.linalg.inv(np.linalg.cholesky(matrix[i]).transpose()),np.linalg.inv(np.linalg.cholesky(matrix[i]))))
#new = np.reshape(invert,(100,4)).transpose()

invert = np.dot(np.linalg.inv(np.linalg.cholesky(matrix).transpose()),np.linalg.inv(np.linalg.cholesky(matrix)))
print invert,np.shape(invert)
"""
plt.plot(masses,np.sqrt(new[0]),color = 'r',marker = '_',label = r'$\omega_b$')
plt.plot(masses,np.sqrt(new[3]),color ='b',marker = '_',label= r'$n_s$')
plt.fill_between(masses,np.sqrt(new[0]),alpha = 0.4,color = 'r')
plt.fill_between(masses,np.sqrt(new[3]),alpha = 0.4,color = 'b')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Mass (eV)')
plt.ylabel(r'$1\sigma$')
plt.title('Forecasted Constraints')
plt.legend(frameon = False, loc = 'upper right')
plt.savefig('test_deriv.png')
subprocess.call('open test_deriv.png',shell = True)
"""
