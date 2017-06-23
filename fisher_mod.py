from numpy import *
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os
import csv
import math
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

#take fresh axicamb and run fixed omega
run_cambAxion4 = False
plot_derivatives = False
plot_wb_derivative = True
plot_wc_derivative = True
plot_fa_derivative = True
plot_ns_derivative = True
plot_h_derivative = True
plot_variance = False
plot_axion_constraints = True
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
#Runs axionCAMB and renames the matterpower output files to reflect the parameters are used for each run
#Each integer (1 through 6) after ./camb params.ini initiates a new parameter value to be set
# 1 = ombh^2 (omegab)     # 4 = mass_axion (ma)
# 2 = omch^2 (omegac)     # 5 = H (H0)
# 3 = omaxh^2 (omegaax)   # 6 = ns (an)
python_directory = os.getcwd()
camb_directory = '/Users/etrott12/Dropbox/emery/axionCAMB'
counter = 0
if run_cambAxion4 == True:
	for j in range(len(masses)):
		for i in range(len(fin_varied[0])):
			os.chdir(camb_directory)
			for k in range(len(fin_varied[0][i])):
				counter += 1
				print 'Starting mass run ',counter,'/',2*len(masses)
				subprocess.call('./camb params.ini 1 %s 2 0.125 3 0.0000001 4 %s 5 67.7 6 0.96'%(fin_varied[0][i][k],masses[j]),shell=True)
				subprocess.call('mv test_matterpower.dat mp/wb%s_m%s.dat'%(counter,j),shell=True)
				subprocess.call('./camb params.ini 1 0.02238 2 %s 3 0.0000001 4 %s 5 67.7 6 0.96'%(fin_varied[1][i][k],masses[j]),shell=True)
				subprocess.call('mv test_matterpower.dat mp/wc%s_m%s.dat'%(counter,j),shell=True)
				subprocess.call('./camb params.ini 1 0.02238 2 %s 3 %s 4 %s 5 67.7 6 0.96'%(central_values[1]-fin_varied[2][i][k],fin_varied[2][i][k],masses[j]),shell=True)
				subprocess.call('mv test_matterpower.dat mp/wa%s_m%s.dat'%(counter,j),shell=True)
				subprocess.call('./camb params.ini 1 0.02238 2 0.125 3 0.0000001 4 %s 5 %s 6 0.96'%(masses[j],fin_varied[3][i][k]),shell=True)
				subprocess.call('mv test_matterpower.dat mp/h%s_m%s.dat'%(counter,j),shell=True)
				subprocess.call('./camb params.ini 1 0.02238 2 0.125 3 0.0000001 4 %s 5 67.7 6 %s'%(masses[j],fin_varied[4][i][k]),shell=True)
				subprocess.call('mv test_matterpower.dat mp/ns%s_m%s.dat'%(counter,j),shell=True)
		subprocess.call('./camb params.ini 1 0.02238 2 0.125 3 0.000001 4 %s 5 67.7 6 0.96'  %(masses[j]), shell = True)
		subprocess.call('mv test_matterpower.dat mp/m%s.dat' %(j), shell = True)
		os.chdir(python_directory)

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

def fin_deriv(k,Pk_param_dict,step_size,param_value):
	param_deriv = {key:[(param_dict[key][1][i]-param_dict[key][0][i])/(step_size*param_value) for i in range(len(param_dict[key][0]))] for key in range(len(k))}
	return param_deriv

#Computes the derivatives with respect to each variable
wb_deriv = {key:[(Pk_wb[key][1][i]-Pk_wb[key][0][i])/(fin_epsilon[0]*central_values[0]) for i in range(len(Pk_wb[key][0]))] for key in range(len(k))}
wc_deriv = {key:[(Pk_wc[key][1][i]-Pk_wc[key][0][i])/(fin_epsilon[1]*central_values[1]) for i in range(len(Pk_wc[key][0]))] for key in range(len(k))}
wa_deriv = {key:[(Pk_wa[key][1][i]-Pk_wa[key][0][i])/(fin_epsilon[2]*central_values[2]) for i in range(len(Pk_wa[key][1]))] for key in range(len(k))}
ns_deriv = {key:[(Pk_ns[key][1][i]-Pk_ns[key][0][i])/(fin_epsilon[3]*central_values[3]) for i in range(len(Pk_ns[key][0]))] for key in range(len(k))}
h_deriv = {key:[(Pk_h[key][1][i]-Pk_h[key][0][i])/(fin_epsilon[4]*central_values[4]) for i in range(len(Pk_h[key][0]))] for key in range(len(k))}

wb_deriv1 = fin_deriv(k,Pk_wb,fin_epsilon[0],central_values[0])
print wb_deriv
print wb_deriv1
"""
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
"""
