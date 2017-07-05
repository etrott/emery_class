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

def run_CLASS(param,value):
    cosmo = Class()
    params = {'output':'tCl,pCl,lCl','lensing':'yes',param:value}
    cosmo.set(params)
    cosmo.compute()
    output = cosmo.lensed_cl(2500)
    return output['tt']

param_list = ['Omega_b','Omega_cdm','H0']
#Sets the order of the polynomial fit (valid from 2 to 8)
poly_fit_order = 6
#Sets the number of points to be run and plotted for each k-value for the polynomial derivative case
num_points = 100
#Sets the central values (around which the derivative will be taken) for each parameter
central_values = [0.1193,67.7,0.02238,0.96]
#Sets the percent change in each parameter for the polynomial derivative case
poly_epsilon = 0.2
#Sets the percent changes in each parameter for the finite difference derivative case
fin_epsilon = [0.065,0.07,0.075,0.085]
fin_matrix = []
poly_matrix = []
for i in range(len(central_values)):
    for j in range(len(fin_epsilon)):
        fin_matrix.append([central_values[i]*(1.0-fin_epsilon[j]),central_values[i]*(1.0+fin_epsilon[j])])
    poly_matrix.append(np.linspace(central_values[i]-central_values[i]*poly_epsilon,central_values[i]+central_values[i]*poly_epsilon,num_points))
fin_matrix = np.reshape(fin_matrix,(len(central_values),len(fin_epsilon),2))
poly_matrix = np.reshape(poly_matrix,(len(central_values),num_points))

cv_matrix = []
for i in range(len(central_values)):
    cv_matrix.append(np.repeat(central_values[i],len(fin_epsilon)))
cv_matrix = np.reshape(cv_matrix,(len(central_values),len(fin_epsilon)))

epsilon_matrix = []
for i in range(len(fin_epsilon)):
    epsilon_matrix.append(np.repeat(fin_epsilon[i],len(fin_epsilon)))
epsilon_matrix = np.reshape(epsilon_matrix,(len(fin_epsilon),len(fin_epsilon)))
print epsilon_matrix,np.shape(epsilon_matrix)
"""
print np.shape(fin_matrix)
print np.shape(poly_matrix)

fin_low = []
fin_high = []
poly = []
for i in range(len(central_values)):
    for k in range(len(fin_epsilon)):
        fin_low.append(run_CLASS(param_list[i],fin_matrix[i][k][0]))
        fin_high.append(run_CLASS(param_list[i],fin_matrix[i][k][1]))
    for j in range(num_points):
        poly.append(run_CLASS(param_list[i],poly_matrix[i][k]))
fin_low = np.reshape(fin_low,(len(central_values),len(fin_epsilon)))
fin_high = np.reshape(fin_high,(len(central_values),len(fin_epsilon)))
poly = np.reshape(poly,len(central_values),num_points))

np.save('fin_low',fin_low)
np.save('fin_high',fin_high)
np.save('poly',poly)

fin_low = np.load('fin_low.npy')
fin_high = np.load('fin_high.npy')
poly = np.load('poly.npy')
"""

def fin_deriv(Pk_high,Pk_low,central_value,step_size):
        param_deriv = (Pk_high-Pk_low)/(central_value*step_size)
        return param_deriv

#fin_deriv = fin_deriv(Pk_high,Pk_low,cv_matrix,epsilon_matrix)
#print fin_deriv,np.shape(fin_deriv)

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
