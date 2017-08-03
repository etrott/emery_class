from classy import Class
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import matplotlib as mpl
from numpy import linalg
from matplotlib.patches import Ellipse
import itertools

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

def run_fisher(param,cv,delta,mass = [1.0],num = 9,deg = 3,run = False, lmax = 2500, plot_cl = False, plot_deriv = False, plot_poly = False,param_index = 6,poly_lmin = 50,poly_lmax = 60):

    def run_CLASS(dict):
        cosmo = Class()
        cosmo.set(dict)
        cosmo.compute()
        output = cosmo.raw_cl(2500)
        #output = cosmo.lensed_cl(3000)
        return output['ell'],output['tt'],output['te'],output['ee']

    mid = (num-1)/2
    cv_floats = []
    cv_strings = []
    param_floats = []
    param_strings = []
    for i in range(len(cv)):
        if isinstance(cv[i],basestring):
            cv_strings.append(cv[i])
            param_strings.append(param[i])
        else:
            cv_floats.append(cv[i])
            param_floats.append(param[i])

    high = np.reshape([cv_floats]*len(cv_floats),(len(cv_floats),len(cv_floats)))
    low = np.reshape([cv_floats]*len(cv_floats),(len(cv_floats),len(cv_floats)))
    span = []
    for i in range(len(high)):
        high[i][i] = high[i][i]+delta[i]
        low[i][i] = low[i][i]-delta[i]
        span.append(np.linspace(low[i][i],high[i][i],num))
    for i in range(len(span)):
        for j in range(len(span[i])):
            if span[i][j] <= 0.:
                #replacing the negative points with positive ones strengthens the constraint by a factor of 3 over setting them to 0
                span[i] = np.linspace(cv_floats[i],high[i][i],num)
                #span[i][j] = 0.
    runs = np.reshape([cv_floats]*len(cv_floats)*num,(len(cv_floats),num,len(cv_floats)))
    for i in range(len(runs)):
        runs[i][:,i] = span[i]

    def save(fname):
        ell = []
        cl_tt = []
        cl_te = []
        cl_ee = []
        counter = 0.
        for i in range(len(runs)):
            for j in range(len(runs[i])):
                counter += 1.
                print counter/np.float((num*len(cv_floats)))
                dict = {param_floats[k]:runs[i][j][k] for k in range(len(param_floats))}
                string_dict = {param_strings[k]:cv_strings[k] for k in range(len(param_strings))}
                dmeff_dict = {'m_dmeff':1.,'cc_dmeff_op':1, 'cc_dmeff_num':1, 'cc_dmeff_n':0, 'cc_dmeff_qm2':0, 'cc_dmeff_scale': 1, 'omega_cdm': 0.0, 'spin_dmeff':0., 'use_helium_dmeff':'yes', 'use_temperature_dmeff':'yes'}
                dict.update(string_dict)
                dict.update(dmeff_dict)
                ell.append(run_CLASS(dict)[0])
                cl_tt.append(run_CLASS(dict)[1])
                cl_te.append(run_CLASS(dict)[2])
                cl_ee.append(run_CLASS(dict)[3])
        ell = np.reshape(ell,(len(cv_floats),num,lmax+1))
        cl_tt = np.reshape(cl_tt,(len(cv_floats),num,lmax+1))
        cl_te = np.reshape(cl_te,(len(cv_floats),num,lmax+1))
        cl_ee = np.reshape(cl_ee,(len(cv_floats),num,lmax+1))
        np.save('data/ell_%s'%(fname),ell)
        np.save('data/cl_tt_%s'%(fname),cl_tt)
        np.save('data/cl_te_%s'%(fname),cl_te)
        np.save('data/cl_ee_%s'%(fname),cl_ee)

    if run == True:
        save('1GeV')

    def load(fname):
        ell = np.load('data/ell_%s.npy'%(fname))
        cl_tt = np.load('data/cl_tt_%s.npy'%(fname))
        cl_te = np.load('data/cl_te_%s.npy'%(fname))
        cl_ee = np.load('data/cl_ee_%s.npy'%(fname))
        return ell,cl_tt,cl_te,cl_ee
    
    ell,cl_tt,cl_te,cl_ee = load('1GeV')

    scale = np.asarray((2.7255**2)*(10**12)).reshape(1,1,1)

    #scaling by temperature
    ell = np.transpose(ell,(0,2,1))
    cl_tt = np.transpose(cl_tt,(0,2,1))
    cl_te = np.transpose(cl_te,(0,2,1))
    cl_ee = np.transpose(cl_ee,(0,2,1))
    cl_tt = cl_tt*scale
    cl_te = cl_te*scale
    cl_ee = cl_ee*scale

    #fits a polynomial to every l for each param
    def poly_fit(x,y,deg):
        poly = []
        for i in range(len(y)):
            for j in range(len(y[i])):
                poly.append(np.polyfit(x[i],y[i][j],deg))
        poly = np.reshape(poly,(len(y),len(y[0]),deg+1))
        return poly

    def plot_polynomial(x,y,deg,param_index,poly_lmin,poly_lmax):
        poly = poly_fit(span,cl_tt,deg)
        test = []
        for i in range(len(poly)):
            for j in range(len(poly[i])):
                for k in range(len(span[i])):
                    test.append(np.poly1d(poly[i][j])(x[i][k]))
        test = np.reshape(test,(len(cv_floats),lmax+1,num))
        for i in range(poly_lmin,poly_lmax):
            fig,ax = plt.subplots()
            ax.plot(x[param_index],test[param_index][i])
            ax.plot(x[param_index],y[param_index][i],'ro')
            plt.savefig('plots/param%s_l%s.pdf' %(param_index,i))
            subprocess.call('open plots/param%s_l%s.pdf' %(param_index,i),shell = True)
        
    if plot_poly == True:
        plot_polynomial(span,cl_tt,deg,param_index,poly_lmin,poly_lmax)
    
    def deriv(cl):
        #computes the polynomial coefficients
        poly_coeff = poly_fit(span,cl,deg)
        #computes the derivative of the polynomial fit
        deriv_coeff = []
        for i in range(len(poly_coeff)):
            for j in range(len(poly_coeff[i])):
                deriv_coeff.append(np.poly1d(np.polyder(poly_coeff[i][j])))
        deriv_coeff = np.reshape(deriv_coeff,(len(poly_coeff),len(poly_coeff[0]),1))
        #evalutes the derivative
        deriv = []
        for i in range(len(deriv_coeff)):
            for j in range(len(deriv_coeff[i])):
                #make it span[i][mid+1] and see if it's the same
                #deriv.append(deriv_coeff[i][j][0](span[i][mid]))
                #print span[i][mid],cv_floats[i]
                #deriv.append(deriv_coeff[i][j][0](cv_floats[i]+span[i][mid+1]))
                deriv.append(deriv_coeff[i][j][0](cv_floats[i]))
        deriv = np.reshape(deriv,(len(deriv_coeff),len(deriv_coeff[0])))
        return deriv

    deriv_tt,deriv_te,deriv_ee = deriv(cl_tt),deriv(cl_te),deriv(cl_ee)

    def plot(x,y,y2,cl_type = 'TT',plot_cl = False, plot_deriv = False):
        y_mod = x*(x+1.0)*y/(2*np.pi)
        y2_mod = x*(x+1.0)*y2/(2*np.pi)
        param_labels = ['\omega_b','\omega_{dmeff}','H_0','n_s','A_s','\\tau','p_{cc}']
        if plot_cl == True:
            fig,ax = plt.subplots()
            ax.plot(x[0][2:],y_mod[0][2:])
            ax.plot(x[0][2:],y2_mod[0][2:])
            plt.xlabel(r'$\ell$')
            plt.ylabel(r'$\ell(\ell+1)C_{\ell}^{%s}/2\pi$ $[\mu \rm{K}^2]$' %(cl_type))
            plt.yscale('log')
            #plt.xscale('log')
            plt.title(r'CMB Power Spectrum')
            plt.savefig('cl.pdf')
            subprocess.call('open cl.pdf',shell = True)
        if plot_deriv == True:
            fig,ax = plt.subplots()
            for i in range(len(cv_floats)):
                #ax.plot(x[i][2:],y_mod[i][2:]*cv_floats[i],label = r'$%s$' %(param_labels[i]))
                if i == 6:
                    ax.plot(x[i][2:],y_mod[i][2:],label = r'$%s$' %(param_labels[i]))
            plt.xlabel(r'$\ell$')
            plt.ylabel(r'$(\ell(\ell+1)/2\pi) \partial{C_{\ell}^{%s}}/\partial{log(p_i)}$ $[\mu \rm{K}^2]$' %(cl_type))
            plt.title(r'CMB Power Spectrum Derivatives')
            plt.legend(frameon = False, loc = 'lower right')
            plt.savefig('plots/deriv_planck.pdf')
            subprocess.call('open plots/deriv_planck.pdf',shell = True)

    if plot_cl == True:
        plot(ell[:,:,mid],cl_tt[:,:,mid],cl_type = 'TT',plot_cl = True)
    if plot_deriv == True:
        plot(ell[:,:,mid],deriv_te,cl_type = 'TE',plot_deriv = True)

    def fisher(deriv,ell,cl,s,theta,fsky,polarization = True):
        ell = ell[:,:,mid]
        cl = cl[:,:,mid]
        s = np.asarray(s).reshape(1,1)
        theta = np.asarray(theta).reshape(1,1)
        n = np.square(s)*np.exp(ell*(ell+1.0)*np.square(theta)/(8.0*np.log(2.0)))
        plot(ell,cl,n,cl_type = 'TT',plot_cl = True)
        #compute cl and derivative matrices
        cl_mat = []
        dmat = []
        if polarization == True:
            cl_mat.append([cl_tt[0,:,mid][2:]+n[0][2:],cl_te[0,:,mid][2:]])
            cl_mat.append([cl_te[0,:,mid][2:],cl_ee[0,:,mid][2:]+(2*n[0][2:])])
            inv = np.linalg.inv(np.reshape(cl_mat,(2,2,lmax-1)).transpose(2,0,1))
            dmat.append([deriv_tt[:,2:],deriv_te[:,2:]])
            dmat.append([deriv_te[:,2:],deriv_ee[:,2:]])
            dmat = np.reshape(dmat,(2,2,len(cv_floats),lmax-1)).transpose(2,3,0,1)
        if polarization == False:
            cl_mat.append([cl_tt[0,:,mid][2:]+n[0][2:]])
            inv = np.linalg.inv(np.reshape(cl_mat,(1,1,lmax-1)).transpose(2,0,1))
            dmat.append([deriv_tt[:,2:]])
            dmat = np.reshape(dmat,(1,1,len(cv_floats),lmax-1)).transpose(2,3,0,1)
        #compute fisher matrix
        fish = []
        for i in range(len(cv_floats)):
            for j in range(len(cv_floats)):
                for k in range(len(ell[0,2:])):
                    fish.append(((2*ell[0][k]+1)*fsky/2)*np.trace(inv[k].dot(dmat[i][k].dot(inv[k]).dot(dmat[j][k]))))
        fish = np.reshape(fish,(len(cv_floats),len(cv_floats),len(ell[0][2:])))
        fish = np.sum(fish,axis = 2)
        return fish

    #planck -- 
    return fisher(deriv,ell,cl_tt,40*np.pi/10800,7*np.pi/10800,0.65,polarization = True)
    #s4 -- return fisher(deriv,ell,cl_tt,1*np.pi/10800,3*np.pi/10800,0.6,polarization = True)
    #cv-limited -- return fisher(deriv,ell,cl_tt,0,3*np.pi/10800,0.50,polarization = True)

#params from table 1 of Planck 2015 results paper
#was running with cc_dmeff_p step size of 5e8
#fisher_matrix = run_fisher(['output','lensing','omega_b','omega_cdm','H0','n_s','A_s','tau_reio'],['tCl,pCl,lCl','yes',0.02222,0.1199,67.26,0.9652,2.199e-9,0.078],[0.00023,0.0022,0.98,0.0062,0.016e-9,0.019])
fisher_matrix = run_fisher(['output','lensing','omega_b','omega_dmeff','H0','n_s','A_s','tau_reio','cc_dmeff_p'],['tCl,pCl','no',0.0222,0.1197,67.31,0.9655,2.2e-9,0.06,0.],[0.00022,0.0012,0.67,0.0097,0.022e-9,0.0006,5e8],run = False)
#print 'Fisher Matrix : ',fisher_matrix

def output_constraints(fisher_matrix,param):
    cov = np.linalg.inv(fisher_matrix)
    for i in range(len(param)):
        constraint = np.sqrt(cov[i][i])
        print '%s : '%(param[i]), constraint

output_constraints(fisher_matrix,['omega_b','omega_dmeff','H0','n_s','A_s','tau_reio','cc_dmeff_p'])

cov = np.linalg.inv(fisher_matrix)
proton_mass = 0.938272 #GeV
dmeff_mass = 0.001 #GeV
mu = (dmeff_mass*proton_mass)/(dmeff_mass+proton_mass)
v = 246 #GeV
sigma = ((cov[6][6]*np.square(mu))/(np.pi*np.power(v,4.)))/np.square(2.5*10**13)
print sigma

"""
cv = [0.0222,0.1197,67.31,0.9655,2.2e-9,0.06,0.]
delta = [0.00022,0.0012,0.67,0.0097,0.022e-9,0.0006,5e8]
alpha_1 = 1.52
alpha_2 = 2.48
alpha_3 = 3.44
hold = [0,1,2,3,4,5,6]
comb = list(itertools.combinations(hold,2))
a2 = []
b2 = []
tan2th = []
for i in range(len(hold)):
    a2.append(((cov[comb[i][0]][comb[i][0]]+cov[comb[i][1]][comb[i][1]])/2.)+np.sqrt((np.square(cov[comb[i][0]][comb[i][0]]-cov[comb[i][1]][comb[i][1]])/4.)+np.square(cov[comb[i][0]][comb[i][1]])))
    b2.append(((cov[comb[i][0]][comb[i][0]]+cov[comb[i][1]][comb[i][1]])/2.)-np.sqrt((np.square(cov[comb[i][0]][comb[i][0]]-cov[comb[i][1]][comb[i][1]])/4.)+np.square(cov[comb[i][0]][comb[i][1]])))
    tan2th.append(2.*cov[comb[i][0]][comb[i][1]]/(cov[comb[i][0]][comb[i][0]]-cov[comb[i][1]][comb[i][1]]))
a = np.sqrt(a2)*alpha_1
b = np.sqrt(b2)*alpha_1
th = np.arctan(tan2th)/2.

def plot_ellipse(param1, param2):
    param_labels = ['\omega_b','\omega_{cdm}','H_0','n_s','A_s','\\tau','p_{cc}']
    plt.figure()
    ax = plt.gca()
    ellipse = Ellipse(xy = (cv[param1],cv[param2]), width = a[param1], height = b[param2], angle = np.rad2deg(th[param1]), alpha = 0.4)
    ax.add_patch(ellipse)
    #plt.ylim(0.,0.0004)
    #plt.xlim(0.0,20000.)
    #plt.xlim(cv[param1]-(alpha_1*delta[param1]),cv[param1]+(alpha_1*delta[param1]))
    #plt.ylim(cv[param2]-(alpha_1*delta[param2]),cv[param2]+(alpha_1*delta[param2]))
    plt.xlabel(r'$%s$' %(param_labels[param1]))
    plt.ylabel(r'$%s$' %(param_labels[param2]))
    plt.savefig('plots/ellipse_%sv%s.pdf' %(param1,param2))
    subprocess.call('open plots/ellipse_%sv%s.pdf' %(param1,param2),shell = True)

#plot_ellipse(5,6)
"""
