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

param_list = ['Omega_b','Omega_cdm','H0']
central_values = [0.022032,0.12038,67.556]
#Fine tuned step sizes to ensure well-behaved derivatives
fin_epsilon_definite = [0.00125,0.005,0.08]
#Converts the step sizes to fractional step sizes
fin_epsilon = [fin_epsilon_definite[i]/central_values[i] for i in range(len(fin_epsilon_definite))]
#Creates a list of masses even space in log-space
masses = np.logspace(-24,-30,10)
ell = np.arange(0,2501).reshape(1,2501)
param_matrix = np.asarray(central_values).reshape(len(param_list),1)
epsilon_matrix = np.asarray(fin_epsilon).reshape(len(param_list),1)
param_values = []
for i in range(len(param_list)):
    param_values.append((central_values[i]*(1.0-fin_epsilon[i]),central_values[i]*(1.0+fin_epsilon[i])))
param_values = np.reshape(param_values,(len(param_list),2))
"""
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
cl_low = np.load('cl_low.npy')
cl_high = np.load('cl_high.npy')

def fin_deriv(cl_high,cl_low,central_value,step_size):
        param_deriv = (cl_high-cl_low)/(central_value*step_size)
        return param_deriv

deriv = fin_deriv(cl_high,cl_low,param_matrix,epsilon_matrix)

param_label = ['\omega_b','\omega_{cdm}','H_0']
fig,ax = plt.subplots()
for i in range(len(param_list)):
    ax.plot(ell[0],deriv[i],label = r'$%s$'%(param_label[i]))
plt.xscale('log')
plt.legend(frameon = False, loc = 'lower right')
plt.savefig('test.png')
subprocess.call('open test.png',shell = True)
