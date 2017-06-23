from numpy import *
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os
import csv
import math
import matplotlib as mpl
from matplotlib.patches import Ellipse

#increase number of mass runs
#(not necessary) plot inverted diagonal element variance (produce same axion variance plot but with parameters fixed instead of marginalized)
#probably physical (see Dan's paper figure 12):
#bump in power spectrum corresponds to (2pi/k) size of horizon at matter-radiation equality - later modes spend a larger amount of their evolution in faster evolving states
#k_eq goes as sqrt(m*H) look up in papers and see if m approx 10^-27

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

#photoz-colors in 5 bands - hard b/c have to deal with instrumentation and error
#as theory to produce a 2d power spectrum - angular positions of galaxies 

#take fresh axicamb and run fixed omega
run_cambAxion4 = False
plot_derivatives = False
plot_wb_derivative = True
plot_wc_derivative = True
plot_fa_derivative = True
plot_ns_derivative = True
plot_h_derivative = True
plot_ellipses = False
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
"""
for j in range(len(masses)):
	for i in range(len(fin_varied[0])):
		os.chdir(camb_directory)
		for k in range(len(fin_varied[0][i])):
			counter += 1
			if counter == 54 and j == 26:
				subprocess.call('./camb params.ini 1 0.02238 2 %s 3 %s 4 %s 5 67.7 6 0.96'%(central_values[1]-fin_varied[2][i][k],fin_varied[2][i][k],masses[j]),shell=True)
				subprocess.call('mv test_matterpower.dat mp/wa%s_m%s.dat'%(counter,j),shell=True)
				print 'omegach^2 = ',central_values[1]-fin_varied[2][i][k]
				print 'omegaaxh^2 = ',fin_varied[2][i][k]
				print 'mass = ',masses[j]
"""
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

#Computes the derivatives with respect to each variable
wb_deriv = {key:[(Pk_wb[key][1][i]-Pk_wb[key][0][i])/(fin_epsilon[0]*central_values[0]) for i in range(len(Pk_wb[key][0]))] for key in range(len(k))}
wc_deriv = {key:[(Pk_wc[key][1][i]-Pk_wc[key][0][i])/(fin_epsilon[1]*central_values[1]) for i in range(len(Pk_wc[key][0]))] for key in range(len(k))}
wa_deriv = {key:[(Pk_wa[key][1][i]-Pk_wa[key][0][i])/(fin_epsilon[2]*central_values[2]) for i in range(len(Pk_wa[key][1]))] for key in range(len(k))}
ns_deriv = {key:[(Pk_ns[key][1][i]-Pk_ns[key][0][i])/(fin_epsilon[3]*central_values[3]) for i in range(len(Pk_ns[key][0]))] for key in range(len(k))}
h_deriv = {key:[(Pk_h[key][1][i]-Pk_h[key][0][i])/(fin_epsilon[4]*central_values[4]) for i in range(len(Pk_h[key][0]))] for key in range(len(k))}

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


"""
#Testing the difference between fisher matrices at different masses
#Small differences between non-wa entries
#Large differences between wa entries
print 'fisher matrix for m = 10^-24'
for m in range(len(matrix_dict[0])):
	if m == 0:
		print 'wb column'      
	if m == 1:
		print 'wc column'
       	if m == 2:
       		print 'wa column'
       	if m == 3:
       		print 'h column'
       	if m == 4:
       		print 'ns column'
	for n in range(len(matrix_dict[0][m])):
		print matrix_dict[0][m][n]
print 'fisher matrix for m = 10^-30'
for m in range(len(matrix_dict[9])):
	if m == 0:
		print 'wb column'
        if m == 1:
                print 'wc column'
        if m == 2:
                print 'wa column'
        if m == 3:
                print 'h column'
        if m == 4:
                print 'ns column'
	for n in range(len(matrix_dict[9][m])):
		print matrix_dict[9][m][n]
"""

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
sub_dict = {key:np.zeros((2,2)) for key in k_wb}
for key in in_range_dict:
	sub_dict[key][0][0] = matrix_dict[key][0][0]
	sub_dict[key][0][1] = matrix_dict[key][0][1]
	sub_dict[key][1][0] = matrix_dict[key][1][0]
	sub_dict[key][1][1] = matrix_dict[key][1][1]
inverted_sub_dict = {key:np.dot(np.linalg.inv(np.linalg.cholesky(sub_dict[key]).transpose()),np.linalg.inv(np.linalg.cholesky(sub_dict[key]))) for key in k_wb}
variance_sub_dict = {i+1:[np.sqrt(inverted_sub_dict[key][i][i]) for key in k_wb] for i in range(0,2)}
"""
"""
#solves the ellipse equation for y given x, a, and b
def ellipse(x,a,b):
        y = np.sqrt(np.square(b)*(1-(np.square(x)/np.square(a))))
        return y
def a_squared(sigma_x,sigma_y,sigma_xy):
        a_squared = ((np.square(sigma_x)+np.square(sigma_y))/2)+np.sqrt((np.square((np.square(sigma_x)-np.square(sigma_y)))/4)+np.square(sigma_xy))
        return a_squared
def b_squared(sigma_x,sigma_y,sigma_xy):
        b_squared = ((np.square(sigma_x)+np.square(sigma_y))/2)-np.sqrt((np.square((np.square(sigma_x)-np.square(sigma_y)))/4)+np.square(sigma_xy))
        return b_squared
def theta(sigma_x,sigma_y,sigma_xy):
        theta = 180*(np.arctan((2*sigma_xy)/(np.square(sigma_x)-np.square(sigma_y))))/(math.pi*2)
        return theta

#computes the degeneracy between each set of parameters 
wb_wc = [[a_squared(variance_dict[1][key],variance_dict[2][key],inverted_dict[key][1][0]),b_squared(variance_dict[1][key],variance_dict[2][key],inverted_dict[key][1][0]),theta(variance_dict[1][key],variance_dict[2][key],inverted_dict[key][1][0])] for key in k_wb]
wb_wa = [[a_squared(variance_dict[1][key-1],variance_dict[3][key-1],inverted_dict[key][2][0]),b_squared(variance_dict[1][key-1],variance_dict[3][key-1],inverted_dict[key][2][0]),theta(variance_dict[1][key-1],variance_dict[3][key-1],inverted_dict[key][2][0])] for key in k_wb]
wb_ns = [[a_squared(variance_dict[1][key-1],variance_dict[4][key-1],inverted_dict[key][3][0]),b_squared(variance_dict[1][key-1],variance_dict[4][key-1],inverted_dict[key][3][0]),theta(variance_dict[1][key-1],variance_dict[4][key-1],inverted_dict[key][3][0])] for key in k_wb]
wb_h = [[a_squared(variance_dict[1][key-1],variance_dict[5][key-1],inverted_dict[key][4][0]),b_squared(variance_dict[1][key-1],variance_dict[5][key-1],inverted_dict[key][4][0]),theta(variance_dict[1][key-1],variance_dict[5][key-1],inverted_dict[key][4][0])] for key in k_wb]
wc_wa = [[a_squared(variance_dict[2][key-1],variance_dict[3][key-1],inverted_dict[key][2][1]),b_squared(variance_dict[2][key-1],variance_dict[3][key-1],inverted_dict[key][2][1]),theta(variance_dict[2][key-1],variance_dict[3][key-1],inverted_dict[key][2][1])] for key in k_wb]
wc_ns = [[a_squared(variance_dict[2][key-1],variance_dict[4][key-1],inverted_dict[key][3][1]),b_squared(variance_dict[2][key-1],variance_dict[4][key-1],inverted_dict[key][3][1]),theta(variance_dict[2][key-1],variance_dict[4][key-1],inverted_dict[key][3][1])] for key in k_wb]
wc_h = [[a_squared(variance_dict[2][key-1],variance_dict[5][key-1],inverted_dict[key][4][1]),b_squared(variance_dict[2][key-1],variance_dict[5][key-1],inverted_dict[key][4][1]),theta(variance_dict[2][key-1],variance_dict[5][key-1],inverted_dict[key][4][1])] for key in k_wb]
wa_ns = [[a_squared(variance_dict[3][key-1],variance_dict[4][key-1],inverted_dict[key][3][2]),b_squared(variance_dict[3][key-1],variance_dict[4][key-1],inverted_dict[key][3][2]),theta(variance_dict[3][key-1],variance_dict[4][key-1],inverted_dict[key][3][2])] for key in k_wb]
wa_h = [[a_squared(variance_dict[3][key-1],variance_dict[5][key-1],inverted_dict[key][4][2]),b_squared(variance_dict[3][key-1],variance_dict[5][key-1],inverted_dict[key][4][2]),theta(variance_dict[3][key-1],variance_dict[5][key-1],inverted_dict[key][4][2])] for key in k_wb]
ns_h = [[a_squared(variance_dict[4][key-1],variance_dict[5][key-1],inverted_dict[key][4][3]),b_squared(variance_dict[4][key-1],variance_dict[5][key-1],inverted_dict[key][4][3]),theta(variance_dict[4][key-1],variance_dict[5][key-1],inverted_dict[key][4][3])] for key in k_wb]

#need to modify to account for the index swap of h and ns in central_values
alpha = 1.52
if plot_ellipses == True:
        images_directory = '/Users/etrott12/Dropbox/emery/axionCAMB/images'
        os.chdir(images_directory)
	fig,ax = plt.subplots()
	for i in range(len(masses)):
		ellipse = mpl.patches.Ellipse(xy = (central_values[0],central_values[1]),width = alpha*wb_wc[i][0],height = alpha*wb_wc[i][1],angle = wb_wc[i][2],fc = 'None',lw = 2,label = masses[i])
		ax.add_patch(ellipse)
	ax.set_xlim(central_values[0]-(alpha*wb_variance[0]),central_values[0]+(alpha*wb_variance[0]))
	ax.set_ylim(central_values[1]-(alpha*wc_variance[0]),central_values[1]+(alpha*wc_variance[0]))
	plt.xlabel('wb')
	plt.ylabel('wc')
	plt.title('Baryon Fraction and CDM Fraction Degeneracies at Multiple Axion Masses')
	plt.savefig('ellipse_wb_wc.png')
	subprocess.call('open ellipse_wb_wc.png',shell = True)
	plt.close('all')

	fig,ax = plt.subplots()
	for i in range(len(masses)):
		ellipse = mpl.patches.Ellipse(xy = (central_values[0],central_values[2]),width = alpha*wb_wa[i][0],height = alpha*wb_wa[i][1],angle = wb_wa[i][2],fc = 'None',lw = 2,label = masses[i])
		ax.add_patch(ellipse)
	ax.set_xlim(central_values[0]-(alpha*wb_variance[0]),central_values[0]+(alpha*wb_variance[0]))
	ax.set_ylim(central_values[2]-(alpha*wa_variance[0]),central_values[2]+(alpha*wa_variance[0]))
	plt.xlabel('wb')
	plt.ylabel('wa')
	plt.title('Baryon Fraction and Axion Fraction Degeneracies at Multiple Axion Masses')
	plt.savefig('ellipse_wb_wa.png')
	subprocess.call('open ellipse_wb_wa.png',shell = True)
	plt.close('all')

	fig,ax = plt.subplots()
	for i in range(len(masses)):
		ellipse = mpl.patches.Ellipse(xy = (central_values[0],central_values[3]),width = alpha*wb_ns[i][0],height = alpha*wb_ns[i][1],angle = wb_ns[i][2],fc = 'None',lw = 2,label = masses[i])
		ax.add_patch(ellipse)
	ax.set_xlim(central_values[0]-(alpha*wb_variance[0]),central_values[0]+(alpha*wb_variance[0]))
	ax.set_ylim(central_values[3]-(alpha*ns_variance[0]),central_values[3]+(alpha*ns_variance[0]))
	plt.xlabel('wb')
	plt.ylabel('ns')
	plt.title('Baryon Fraction and Ns Degeneracies at Multiple Axion Masses')
	plt.savefig('ellipse_wb_ns.png')
	subprocess.call('open ellipse_wb_ns.png',shell = True)
	plt.close('all')

	fig,ax = plt.subplots()
	for i in range(len(masses)):
		ellipse = mpl.patches.Ellipse(xy = (central_values[0],central_values[4]),width = alpha*wb_h[i][0],height = alpha*wb_h[i][1],angle = wb_h[i][2],fc = 'None',lw = 2,label = masses[i])
		ax.add_patch(ellipse)
	ax.set_xlim(central_values[0]-(alpha*wb_variance[0]),central_values[0]+(alpha*wb_variance[0]))
	ax.set_ylim(central_values[4]-(alpha*h_variance[0]),central_values[4]+(alpha*h_variance[0]))
	plt.xlabel('wb')
	plt.ylabel('h')
	plt.title('Baryon Fraction and Hubble Constant Degeneracies at Multiple Axion Masses')
	plt.savefig('ellipse_wb_h.png')
	subprocess.call('open ellipse_wb_h.png',shell = True)
	plt.close('all')

	fig,ax = plt.subplots()
	for i in range(len(masses)):
		ellipse = mpl.patches.Ellipse(xy = (central_values[1],central_values[2]),width = alpha*wc_wa[i][0],height = alpha*wc_wa[i][1],angle = wc_wa[i][2],fc = 'None',lw = 2,label = masses[i])
		ax.add_patch(ellipse)
	ax.set_xlim(central_values[1]-(alpha*wc_variance[0]),central_values[1]+(alpha*wc_variance[0]))
	ax.set_ylim(central_values[2]-(alpha*wa_variance[0]),central_values[2]+(alpha*wa_variance[0]))
	plt.xlabel('wc')
	plt.ylabel('wa')
	plt.title('CDM Fraction and Axion Fraction Degeneracies at Multiple Axion Masses')
	plt.savefig('ellipse_wc_wa.png')
	subprocess.call('open ellipse_wc_wa.png',shell = True)
	plt.close('all')

	fig,ax = plt.subplots()
	for i in range(len(masses)):
		ellipse = mpl.patches.Ellipse(xy = (central_values[1],central_values[3]),width = alpha*wc_ns[i][0],height = alpha*wc_ns[i][1],angle = wc_ns[i][2],fc = 'None',lw = 2,label = masses[i])
		ax.add_patch(ellipse)
	ax.set_xlim(central_values[1]-(alpha*wc_variance[0]),central_values[1]+(alpha*wc_variance[0]))
	ax.set_ylim(central_values[3]-(alpha*ns_variance[0]),central_values[3]+(alpha*ns_variance[0]))
	plt.xlabel('wc')
	plt.ylabel('ns')
	plt.title('CDM Fraction and Ns Degeneracies at Multiple Axion Masses')
	plt.savefig('ellipse_wc_ns.png')
	subprocess.call('open ellipse_wc_ns.png',shell = True)
	plt.close('all')

	fig,ax = plt.subplots()
	for i in range(len(masses)):
		ellipse = mpl.patches.Ellipse(xy = (central_values[1],central_values[4]),width = alpha*wc_h[i][0],height = alpha*wc_h[i][1],angle = wc_h[i][2],fc = 'None',lw = 2,label = masses[i])
		ax.add_patch(ellipse)
	ax.set_xlim(central_values[1]-(alpha*wc_variance[0]),central_values[1]+(alpha*wc_variance[0]))
	ax.set_ylim(central_values[4]-(alpha*h_variance[0]),central_values[4]+(alpha*h_variance[0]))
	plt.xlabel('wc')
	plt.ylabel('h')
	plt.title('CDM Fraction and Hubble Constant Degeneracies at Multiple Axion Masses')
	plt.savefig('ellipse_wc_h.png')
	subprocess.call('open ellipse_wc_h.png',shell = True)
	plt.close('all')

	fig,ax = plt.subplots()
	for i in range(len(masses)):
		ellipse = mpl.patches.Ellipse(xy = (central_values[2],central_values[3]),width = alpha*wa_ns[i][0],height = alpha*wa_ns[i][1],angle = wa_ns[i][2],fc = 'None',lw = 2,label = masses[i])
		ax.add_patch(ellipse)
	ax.set_xlim(central_values[2]-(alpha*wa_variance[0]),central_values[2]+(alpha*wa_variance[0]))
	ax.set_ylim(central_values[3]-(alpha*ns_variance[0]),central_values[3]+(alpha*ns_variance[0]))
	plt.xlabel('wa')
	plt.ylabel('ns')
	plt.title('Axion Fraction and Ns Degeneracies at Multiple Axion Masses')
	plt.savefig('ellipse_wa_ns.png')
	subprocess.call('open ellipse_wa_ns.png',shell = True)
	plt.close('all')

	fig,ax = plt.subplots()
	for i in range(len(masses)):
		ellipse = mpl.patches.Ellipse(xy = (central_values[2],central_values[4]),width = alpha*wa_h[i][0],height = alpha*wa_h[i][1],angle = wa_h[i][2],fc = 'None',lw = 2,label = masses[i])
		ax.add_patch(ellipse)
	ax.set_xlim(central_values[2]-(wa_variance[0]*alpha),central_values[2]+(wa_variance[0]*alpha))
	ax.set_ylim(central_values[4]-(h_variance[0]*alpha),central_values[4]+(h_variance[0]*alpha))
	plt.xlabel('wa')
	plt.ylabel('h')
	plt.title('Axion Fraction  and Hubble Constant Degeneracies at Multiple Axion Masses')
	plt.savefig('ellipse_wa_h.png')
	subprocess.call('open ellipse_wa_h.png',shell = True)
	plt.close('all')

	fig,ax = plt.subplots()
	for i in range(len(masses)):
		ellipse = mpl.patches.Ellipse(xy = (central_values[3],central_values[4]),width = alpha*ns_h[i][0],height = alpha*ns_h[i][1],angle = ns_h[i][2],fc = 'None',lw = 2,label = masses[i])
		ax.add_patch(ellipse)
	ax.set_xlim(central_values[3]-(ns_variance[0]*alpha),central_values[3]+(ns_variance[0]*alpha))
	ax.set_ylim(central_values[4]-(h_variance[0]*alpha),central_values[4]+(h_variance[0]*alpha))
	plt.xlabel('ns')
	plt.ylabel('h')
	plt.title('Ns  and Hubble Constant Degeneracies at Multiple Axion Masses')
	plt.savefig('ellipse_ns_h.png')
	subprocess.call('open ellipse_ns_h.png',shell = True)
	plt.close('all')
"""
dan_mass = [9.8e-31,1.76e-30,4.11e-30,7.23e-29,2.21e-27,7.14e-27,1.2e-26,1.92e-26,2.96e-26,4.71e-26,6.2e-26,1.14e-25,1.33e-25,1.59e-25,2.18e-25,1.0e-24]
dan_sigma = [6.53E-02,6.75e-02,6.98e-02,6.86e-02,6.96e-02,8.01e-02,1.15e-01,2.17e-01,3.54e-01,4.99e-01,6.09e-01,6.63e-01,7.46e-01,8.16e-01,8.87e-01,1.0]
if plot_axion_constraints == True:
#        images_directory = '/Users/etrott12/Dropbox/emery/axionCAMB/images'
#        os.chdir(images_directory)
	fig,ax = plt.subplots()
	ax.fill_between(masses,wa_variance,0,color = 'b',alpha = 0.5,label = r'Forecasted SDSS Constraint')
        ax.plot(masses,wa_variance,color = 'b',label = r'Forecasted SDSS Constraint')
	ax.fill_between(dan_mass,dan_sigma,0,color = 'r',alpha = 0.3,label = r'CMB constraint')
	ax.plot(dan_mass,dan_sigma,color = 'r',label = r'CMB constraint')
	plt.xlim(10e-31,10e-25)
	plt.ylim(10e-4,10e0)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$m_a$ (eV)')
	plt.ylabel(r'$f_a$ $1\sigma$ Variance')
	plt.title(r'Forecasted SDDS $f_a$ Constraints')
	plt.legend(frameon = False, loc = 'upper left')
	plt.savefig('plots/axion_variance.pdf')
	subprocess.call('open plots/axion_variance.pdf',shell = True)
	plt.close('all')


if plot_variance == True:
#        images_directory = '/Users/etrott12/Dropbox/emery/axionCAMB/images'
#        os.chdir(images_directory)
	fig,ax = plt.subplots()
	ax.plot(masses,wb_variance,color = 'b',label = 'wb')
	ax.plot(masses,wc_variance,color = 'r',label = 'wc')
	ax.plot(masses,wa_variance,color = 'g',label = 'fa')
	ax.plot(masses,ns_variance,color = 'k',label = 'ns')
	ax.plot(masses,h_variance,color = 'c',label = 'H')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Axion Mass (eV)')
	plt.ylabel('One-Sigma Variance')
	plt.title('One-Sigma Parameter Constraints from SDSS z = 0.3')
	plt.legend(frameon = False, loc = 'upper right')
	plt.savefig('plots/variance.png')
	subprocess.call('open plots/variance.png',shell = True)
	plt.close('all')

if plot_derivatives == True:
#	images_directory = '/Users/etrott12/Dropbox/emery/axionCAMB/images'
#	os.chdir(images_directory)
	if plot_wb_derivative == True:
		fig,ax = plt.subplots()
		for key in k_wb:
			ax.plot(k_wb[key][1],wb_deriv[key],label = 'm = %s' %(masses[key]))
			ax.set_xlim(0.0001,1.0)
                        ax.set_xscale('log')
                        plt.legend(frameon = False, loc = 'lower left')
		plt.savefig('plots/deriv_wb.png')
		subprocess.call('open plots/deriv_wb.png',shell = True)
		plt.close('all')
	if plot_wc_derivative == True:
		fig,ax = plt.subplots()
		for key in k_wb:
                        ax.plot(k_wc[key][0],wc_deriv[key],label = 'm = %s' %(masses[key]))
                        ax.set_xlim(0.0001,1.0)
                        ax.set_xscale('log')
                        plt.legend(frameon = False, loc = 'lower left')
		plt.savefig('plots/deriv_wc.png')
		subprocess.call('open plots/deriv_wc.png',shell = True)
		plt.close('all')
	if plot_fa_derivative == True:
		fig,ax = plt.subplots()
		for key in k_wb:
                        ax.plot(k_wa[key][0],wa_deriv[key],label = 'm = %s' %(masses[key]))
                        ax.set_xlim(0.0001,1.0)
                        ax.set_xscale('log')
                        plt.legend(frameon = False, loc = 'lower left')
		plt.savefig('plots/deriv_wa.png')
		subprocess.call('open plots/deriv_wa.png',shell = True)
		plt.close('all')
	if plot_ns_derivative == True:
		fig,ax = plt.subplots()
		for key in k_wb:
                        ax.plot(k_ns[key][0],ns_deriv[key],label = 'm = %s' %(masses[key]))
                        ax.set_xlim(0.0001,1.0)
                        ax.set_xscale('log')
                        plt.legend(frameon = False, loc = 'lower left')
		plt.savefig('plots/deriv_ns.png')
		subprocess.call('open plots/deriv_ns.png',shell = True)
		plt.close('all')
	if plot_h_derivative == True:
		fig,ax = plt.subplots()
		for key in k_wb:
                        ax.plot(k_h[key][0],h_deriv[key],label = 'm = %s' %(masses[key]))
                        ax.set_xlim(0.0001,1.0)
                        ax.set_xscale('log')
                        plt.legend(frameon = False, loc = 'lower left')
		plt.savefig('plots/deriv_h.png')
		subprocess.call('open plots/deriv_h.png',shell = True)
		plt.close('all')
