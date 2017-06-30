from classy import Class
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import subprocess

def run_CLASS(param,value):
    cosmo = Class()
    params = {'output':'tCl,pCl,lCl','lensing':'yes',param:value}
    cosmo.set(params)
    cosmo.compute()
    output = cosmo.lensed_cl(2500)
    return output['tt']

param_list = ['omega_b','omega_cdm','h','n_s']
central_values = [0.02238,0.125,67.7,0.96]
#Fine tuned step sizes to ensure well-behaved derivatives
fin_epsilon_definite = [0.00125,0.0051,0.08,0.025]
#Converts the step sizes to fractional step sizes
fin_epsilon = [fin_epsilon_definite[i]/central_values[i] for i in range(len(fin_epsilon_definite))]
#Creates upper and lower values for each parameter to obtain accurate derivatives
#fin_varied = {key:[] for key in range(len(central_values))}
#for key in fin_varied:
#    fin_varied[key].append([central_values[key]*(1.0-fin_epsilon[key]),central_values[key]*(1.0+fin_epsilon[key])])
#Creates a list of masses even space in log-space
#masses = np.logspace(-24,-30,10)
#print fin_varied

ell = np.arange(0,2501).reshape(1,2501)

param_matrix = np.asarray(central_values).reshape(len(param_list),1)
epsilon_matrix = np.asarray(fin_epsilon).reshape(len(param_list),1)

param_values = []
for i in range(len(param_list)):
    param_values.append((central_values[i]*(1.0-fin_epsilon[i]),central_values[i]*(1.0+fin_epsilon[i])))
param_values = np.reshape(param_values,(len(param_list),2))
print list

cl_low = []
cl_high = []
for i in range(len(param_list)):
    cl_low.append(run_CLASS(param_list[i],param_values[i][0]))
    cl_high.append(run_CLASS(param_list[i],param_values[i][1]))
cl_low = np.reshape(cl_low,(len(param_list),2501))
cl_high = np.reshape(cl_high,(len(param_list),2501))

def fin_deriv(cl_high,cl_low,central_value,step_size):
        param_deriv = (cl_high-cl_low)/(central_value*step_size)
        return param_deriv

deriv = fin_deriv(cl_high,cl_low,param_matrix,epsilon_matrix)

"""
lensed_cl_low = run_CLASS('h',0.62)
lensed_cl_high = run_CLASS('h',0.72)

factor = lensed_cl['ell']*(lensed_cl['ell']+1)/(2*np.pi)
fig,ax = plt.subplots()
ax.loglog(lensed_cl['ell'],factor*lensed_cl['tt'],label = '$h = 0.62$')
ax.loglog(lensed_cl_high['ell'],factor*lensed_cl_high['tt'],label = '$h = 0.72$')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
plt.xlim(2,3000)
plt.ylim(10**-11,10**-9)
plt.legend(frameon = False,loc = 'lower left')
plt.savefig('test.png')
subprocess.call('open test.png',shell = True)
"""
