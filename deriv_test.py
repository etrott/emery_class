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

one_sided_deriv = False
two_sided_deriv = True
#Sets the order of the polynomial fit (valid from 2 to 8)
poly_fit_order = 6
#Sets the number of points to be run and plotted for each k-value for the polynomial derivative case
num_points = 100
#Sets the central values (around which the derivative will be taken) for each parameter
central_values = [0.1193,67.7,0.02238,0.96]
#Creates a dictionary to hold the polynomial derivative parameter values that cambAxion4 will be run at
poly_varied = {key: [] for key in range(1,len(central_values)+1)}
#Sets the percent change in each parameter for the polynomial derivative case
poly_epsilon = 0.4
#Fills the 'poly_varied' dictionary with polynomial derivative parameters values that cambAxion4 will be run at -- num_points values centered around the central parameter value
for key in poly_varied:
	poly_varied[key].append(np.linspace(central_values[key-1]-central_values[key-1]*poly_epsilon,central_values[key-1]+central_values[key-1]*poly_epsilon,num_points))
#Sets the percent changes in each parameter for the finite difference derivative case
fin_epsilon = [0.065,0.07,0.075,0.085]
#Creates a dictionary to hold the finite difference derivative parameter values that cambAxion4 will be run at
fin_varied = {key: [] for key in range(1,len(central_values)+1)}
#Fills the 'fin_varied' dictionary with the finite difference derivative parameter values that cambAxion4 will be run at
for key in fin_varied:
	for i in range(len(fin_epsilon)):
		if one_sided_deriv == True:
			fin_varied[key].append([central_values[key-1],central_values[key-1]*(1.0+fin_epsilon[i])])
		if two_sided_deriv == True:
			fin_varied[key].append([central_values[key-1]*(1.0-fin_epsilon[i]),central_values[key-1]*(1.0+fin_epsilon[i])])

"""
#Runs cambAxion4 and renames the matterpower output files to reflect the parameters are used for each run         
#Each integer (1 through 9) after ./camb emery_params_standard.ini initiates a new parameter value to be set      
# 1 = omch^2 (omegac)     6 = tau (optical_depth)                                                                 
# 2 = H (H0)              7 = ns (an)                                                                             
# 3 = ombh^2 (omegab)     8 = As (ScalarPowerAmp(1))                                                              
# 4 = omaxh^2 (omegaax)   9 = Al (Alens)                                                                          
# 5 = mass_axion (ma) 
python_directory = os.getcwd()
camb_directory = '/Users/etrott12/Dropbox/emery/cambAxion4'
for i in range(num_points):
	os.chdir(camb_directory)
	if run_poly_camb == True:
		print 'Starting cambAxion4 poly_camb run ',i+1,'/',num_points
		subprocess.call('./camb emery_params_poly.ini 1 0.1134 2 %s 3 0.02238 4 10.e-27 5 0.00597 6 0.082 7 0.96 8 2.2132e-9 9 1.00' %(poly_varied[param_number][0][i]), shell = True)
		subprocess.call('mv test_matterpower_xxx.dat poly_h%s.dat' %(i+1), shell = True)
	os.chdir(python_directory)
counter = 0
for i in range(len(fin_varied[1])):
	os.chdir(camb_directory)
	if run_fin_camb == True:
		for k in range(len(fin_varied[param_number][i])):
			counter += 1
			print 'Starting cambAxion4 fin_camb run ',counter,'/',len(fin_varied[1])*2
			subprocess.call('./camb emery_params_poly.ini 1 0.1134 2 %s 3 0.02238 4 10.e-27 5 0.00597 6 0.082 7 0.96 8 2.2132e-9 9 1.00' %(fin_varied[param_number][i][k]), shell = True)
			subprocess.call('mv test_matterpower_xxx.dat fin_h%s.dat' %(counter), shell = True)
	os.chdir(python_directory)
"""

def load_poly(param_string_list,param_index_list):
        k_poly = []
        Pk_poly = []
        for j in range(len(param_string_list)):
                for i in range(len(param_index_list)):
                        k_poly.append(loadtxt('/Users/etrott12/Dropbox/emery/cambAxion4/poly_%s%s.dat' %(param_string_list[j],param_index_list[i]), usecols = [0]))
                        Pk_poly.append(loadtxt('/Users/etrott12/Dropbox/emery/cambAxion4/poly_%s%s.dat' %(param_string_list[j],param_index_list[i]), usecols = [1]))
#        k_poly = np.reshape(k_poly,(len(param_string_list),14,551))
#        Pk_poly = np.reshape(Pk_poly,(len(param_string_list),14,551))
        return k_poly,Pk_poly

def load_fin(param_string_list,param_index_list):
        k_fin = []
        Pk_fin = []
        for j in range(len(param_string_list)):
                for i in range(len(param_index_list)):
                        k_fin.append(loadtxt('/Users/etrott12/Dropbox/emery/cambAxion4/fin_%s%s.dat' %(param_string_list[j],param_index_list[i]), usecols = [0]))
                        Pk_fin.append(loadtxt('/Users/etrott12/Dropbox/emery/cambAxion4/fin_%s%s.dat' %(param_string_list[j],param_index_list[i]), usecols = [1]))
#        k_fin = np.reshape(k_fin,(len(param_string_list),14,551))
#        Pk_fin = np.reshape(Pk_fin,(len(param_string_list),14,551))
        return k_fin,Pk_fin

k_poly,Pk_poly = load_poly(['h'],np.arange(1,15))
k_low,Pk_low = np.asarray(load_fin(['h'],np.arange(1,15,2)))
k_high,Pk_high = np.asarray(load_fin(['h'],np.arange(2,15,2)))

#print Pk_low,Pk_high
#print central_values[1], fin_epsilon[1]

def fin_deriv(Pk_high,Pk_low,central_value,step_size):
        param_deriv = (Pk_high-Pk_low)
	#/(central_value*step_size)
        return param_deriv

fin_deriv = fin_deriv(Pk_high,Pk_low,central_values[1],fin_epsilon[1])
#print fin_deriv,np.shape(fin_deriv)

k_sorted = np.transpose(fin_deriv)
#plt.plot(k_fin,k_sorted)

"""
k_max_list = []
for i in range(len(k_H_fin[1])):
	if k_H_fin[1][i] < 1.1:
		k_max_list.append(i)
last_k = np.amax(k_max_list)
k_range = [k_H_fin[1][i] for i in range(last_k)]

#Creates and fills a dictionary whose keys correspond to the index fo each k-value and whose values correspond to the derivatives of each fin_epsilon pair at the corresponding k-value
fin_deriv_H = {key: [] for key in range(1,len(fin_varied[1])+1)}
for key in fin_deriv_H:
	for i in range(last_k):
		if one_sided_deriv == True:
			fin_deriv_H[key].append((Pk_H_fin[2*key][i]-Pk_H_fin[2*key-1][i])/(fin_epsilon[key-1]))
		if two_sided_deriv == True:
			fin_deriv_H[key].append((Pk_H_fin[2*key][i]-Pk_H_fin[2*key-1][i])/(2*fin_epsilon[key-1]))

#Creates a dictionary whose keys correspond the the index of each k-value and whose values correspond to all of the Pk-values at the k-value index indicated by the key
H_poly_dict = {key: [] for key in range(1,last_k+1)}
for k in range(last_k):
	for i in range(1,num_points+1):
		H_poly_dict[k+1].append(Pk_H_poly[i][k])

#Creates a list that will hold all of the coefficients from the first order term in the polynomial
first_deriv_coeff = []
#Creates a list that will hold all of the 1-r^2 values (one for each k) to check how well the polynomial fit approximates the cambAxion4 output points
one_minus_rsquared_poly = []
#Converts the parameter values to be relative to the central value for polynomial fitting purposes
H_poly_values = [poly_varied[param_number][0][i]-central_values[param_number-1] for i in range(len(poly_varied[param_number][0]))]

for key in H_poly_dict:
	#computes a polynomial fit for each k-value and outputs the coefficients in descending order of powers
        poly_fit = np.polyfit(H_poly_values,H_poly_dict[key],poly_fit_order)
	first_deriv_coeff.append(poly_fit[poly_fit_order-1]*central_values[param_number-1])
	#creates a list of x-values that the polynomial fit will be computed at
	x_poly_list = np.linspace(H_poly_values[0],H_poly_values[len(H_poly_values)-1],num_points)
	#computes the polynomial fit value at each x-value in x_list
        if poly_fit_order == 8:
                y_poly_list = [poly_fit[0]*i**8+poly_fit[1]*i**7+poly_fit[2]*i**6+poly_fit[3]*i**5+poly_fit[4]*i**4+poly_fit[5]*i**3+poly_fit[6]*i**2+poly_fit[7]*i+poly_fit[8] for i in x_poly_list]
        if poly_fit_order == 7:
                y_poly_list = [poly_fit[0]*i**7+poly_fit[1]*i**6+poly_fit[2]*i**5+poly_fit[3]*i**4+poly_fit[4]*i**3+poly_fit[5]*i**2+poly_fit[6]*i+poly_fit[7] for i in x_poly_list]
	if poly_fit_order == 6:
		y_poly_list = [poly_fit[0]*i**6+poly_fit[1]*i**5+poly_fit[2]*i**4+poly_fit[3]*i**3+poly_fit[4]*i**2+poly_fit[5]*i+poly_fit[6] for i in x_poly_list]
	if poly_fit_order == 5:
		y_poly_list = [poly_fit[0]*i**5+poly_fit[1]*i**4+poly_fit[2]*i**3+poly_fit[3]*i**2+poly_fit[4]*i+poly_fit[5] for i in x_poly_list]
	if poly_fit_order == 4:
		y_poly_list = [poly_fit[0]*i**4+poly_fit[1]*i**3+poly_fit[2]*i**2+poly_fit[3]*i+poly_fit[4] for i in x_poly_list]
	if poly_fit_order == 3:
		y_poly_list = [poly_fit[0]*i**3+poly_fit[1]*i**2+poly_fit[2]*i+poly_fit[3] for i in x_poly_list]
	if poly_fit_order == 2:
		y_poly_list = [poly_fit[0]*i**2+poly_fit[1]*i+poly_fit[2] for i in x_poly_list]
	#plots the Pk-values for each k-value as a function of the parameter value and overplots the polynomial fit
	if plot_poly_fits == True:
		fig,ax = plt.subplots()
		print 'Plotting k = ',key,' out of ',len(H_poly_dict)
		for i in range(len(H_poly_dict[key])):
			ax.plot(H_poly_values[i],H_poly_dict[key][i],color = 'k', marker = '.', label = 'k = ',)
		ax.plot(x_poly_list,y_poly_list,color = 'b')
		ax.set_xlabel('H')
		ax.set_ylabel('P(k)')
		ax.set_ylim(np.amin(H_poly_dict[key])-0.1*np.amin(H_poly_dict[key]),np.amax(H_poly_dict[key])+0.1*np.amax(H_poly_dict[key]))
		plt.savefig('/Users/etrott12/Dropbox/emery/CAMB-Feb2015/python/plots/k%s.png' %(key))
		plt.close('all')
	#ybar is the average value Pk for each value of k over the specified parameter value range
	ybar_poly = np.sum(H_poly_dict[key])/len(H_poly_dict)
	#sstot is the sum of the squares of the difference between ybar (the average Pk) and each individual Pk-value
	sstot_poly = 0
	for i in range(len(H_poly_dict[key])):
		sstot_poly += np.square(H_poly_dict[key][i]-ybar_poly)
	#ssres is the sum of the squares of the difference between ybar (the average Pk) and the polynomial fit
	ssres_poly = 0
	for i in range(len(x_poly_list)):
		ssres_poly += np.square(H_poly_dict[key][i]-y_poly_list[i])
	one_minus_rsquared_poly.append(ssres_poly/sstot_poly)
	
ybar_fin_epsilon = np.sum(first_deriv_coeff)/len(first_deriv_coeff)
sstot_fin_epsilon = 0
for i in range(len(first_deriv_coeff)):
	sstot_fin_epsilon += np.square(first_deriv_coeff[i]-ybar_fin_epsilon)
ssres_fin_epsilon = {key: [] for key in range(1,len(fin_epsilon)+1)}
one_minus_rsquared_fin_epsilon = {key: [] for key in range(1,len(fin_epsilon)+1)}
for i in range(1,len(fin_deriv_H)+1):
	for k in range(len(fin_deriv_H[i])):
		ssres_fin_epsilon[i].append(np.square(first_deriv_coeff[k]-fin_deriv_H[i][k]))


sum = []
for key in one_minus_rsquared_fin_epsilon:
	holder = 0
	for i in range(len(ssres_fin_epsilon[key])):
		holder += ssres_fin_epsilon[key][i]/sstot_fin_epsilon
	sum.append(holder)

"""
