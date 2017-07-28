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

def run_fisher(param,cv,delta,mass = [1.0],num = 7,deg = 3,run = False, lmax = 2500):

    def run_CLASS(dict):
        cosmo = Class()
        cosmo.set(dict)
        cosmo.compute()
        output = cosmo.raw_cl(lmax)
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
    runs = np.reshape([cv_floats]*len(cv_floats)*num,(len(cv_floats),num,len(cv_floats)))
    for i in range(len(runs)):
        runs[i][:,i] = span[i]

    def save():
        ell = []
        cl_tt = []
        cl_te = []
        cl_ee = []
        for i in range(len(runs)):
            for j in range(len(runs[i])):
                dict = {param_floats[k]:runs[i][j][k] for k in range(len(param_floats))}
                string_dict = {param_strings[k]:cv_strings[k] for k in range(len(param_strings))}
                dict.update(string_dict)
                ell.append(run_CLASS(dict)[0])
                cl_tt.append(run_CLASS(dict)[1])
                cl_te.append(run_CLASS(dict)[2])
                cl_ee.append(run_CLASS(dict)[3])
        ell = np.reshape(ell,(len(cv_floats),num,lmax+1))
        cl_tt = np.reshape(cl_tt,(len(cv_floats),num,lmax+1))
        cl_te = np.reshape(cl_te,(len(cv_floats),num,lmax+1))
        cl_ee = np.reshape(cl_ee,(len(cv_floats),num,lmax+1))
        np.save('data/ell_mat',ell)
        np.save('data/cl_tt_mat',cl_tt)
        np.save('data/cl_te_mat',cl_te)
        np.save('data/cl_ee_mat',cl_ee)

    if run == True:
        save()

    def load():
        ell = np.load('data/ell_mat.npy')
        cl_tt = np.load('data/cl_tt_mat.npy')
        cl_te = np.load('data/cl_te_mat.npy')
        cl_ee = np.load('data/cl_ee_mat.npy')
        return ell,cl_tt,cl_te,cl_ee
    
    ell,cl_tt,cl_te,cl_ee = load()

    scale = np.asarray((2.7255**2)*(10**12)).reshape(1,1,1)

    ell = np.transpose(ell,(0,2,1))
    cl_tt = np.transpose(cl_tt,(0,2,1))
    cl_te = np.transpose(cl_te,(0,2,1))
    cl_ee = np.transpose(cl_ee,(0,2,1))
    cl_tt = cl_tt*scale
    cl_te = cl_te*scale
    cl_ee = cl_ee*scale

    def poly_fit(x,y,deg):
        poly = []
        for i in range(len(y)):
            for j in range(len(y[i])):
#                poly.append(np.poly1d(np.polyfit(x[i],y[i][j],deg)))
                poly.append(np.polyfit(x[i],y[i][j],deg))
        poly = np.reshape(poly,(len(y),len(y[0]),deg+1))
        return poly
    """
        poly = np.reshape(poly,(len(y),len(y[0])))
        test = []
        for i in range(len(poly)):
            for j in range(len(poly[i])):
                for k in range(len(span[i])):
                    test.append(poly[i][j](span[i][k]))
        test = np.reshape(test,(3,2501,25))
        return test
        return poly
    a = poly_fit(span,cl_tt,deg)
    print a,np.shape(a)
    fig,ax = plt.subplots()
    ax.plot(span[2],a[2][2500],'b')
    ax.plot(span[2],cl_tt[2][2500],'ro')
    plt.show()
    """

    def deriv(cl):
        poly_coeff = poly_fit(span,cl,deg)
        deriv_coeff = []
        for i in range(len(poly_coeff)):
            for j in range(len(poly_coeff[i])):
                deriv_coeff.append(np.poly1d(np.polyder(poly_coeff[i][j])))
        deriv_coeff = np.reshape(deriv_coeff,(len(poly_coeff),len(poly_coeff[0]),1))
        deriv = []
        for i in range(len(deriv_coeff)):
            for j in range(len(deriv_coeff[i])):
                deriv.append(deriv_coeff[i][j][0](span[i][mid]))
        deriv = np.reshape(deriv,(len(deriv_coeff),len(deriv_coeff[0])))
        return deriv

    deriv_tt,deriv_te,deriv_ee = deriv(cl_tt),deriv(cl_te),deriv(cl_ee)

    def plot(x,y,deriv = False):
        y_mod = x*(x+1.0)*y/(2*np.pi)
        param_labels = ['\omega_b','\omega_{cdm}','H_0','n_s','A_s','\\tau']
        fig,ax = plt.subplots()
        if deriv == False:
            ax.plot(x[0][2:],y_mod[0][2:],label = '$TT$')
        if deriv == True:
#            ax.plot(x[2][2:],y_mod[2][2:]*67.8)
            for i in range(len(y_mod)):
                ax.plot(x[0][2:],y_mod[i][2:]*cv_floats[i],label = r'$%s$' %(param_labels[i]))
        plt.xlabel(r'$\ell$')
        if deriv == False:
            plt.ylabel(r'$\ell(\ell+1)C_{\ell}/2\pi$ $[\mu \rm{K}^2]$')
            plt.title(r'CMB Power Spectrum')
            plt.legend(frameon = False, loc = 'upper right')
            plt.savefig('cl.pdf')
            subprocess.call('open cl.pdf',shell = True)
        if deriv == True:
            plt.ylabel(r'$\partial{C_l}/\partial{log(p_i)}$ $[\mu \rm{K}^2]$')
            plt.title(r'CMB Power Spectrum Derivatives')
            plt.legend(frameon = False, loc = 'upper right')
            plt.savefig('plots/deriv_planck.pdf')
            subprocess.call('open plots/deriv_planck.pdf',shell = True)

#    plot(ell[:,:,mid],cl_tt[:,:,mid])
#    plot(ell[:,:,mid],deriv_tt,deriv = True)

    def fisher(deriv,ell,cl,s,theta,fsky):
        ell = ell[:,:,mid]
        cl = cl[:,:,mid]
        s = np.asarray(s).reshape(1,1)
        theta = np.asarray(theta).reshape(1,1)
        n = np.square(s)*np.exp(ell*(ell+1.0)*np.square(theta)/(8.0*np.log10(2.0)))
        cl_mat = []
        cl_mat.append([cl_tt[0,:,mid][2:]+n[0][2:],cl_te[0,:,mid][2:]])
        cl_mat.append([cl_te[0,:,mid][2:],cl_ee[0,:,mid][2:]+(2*n[0][2:])])
        inv = np.linalg.inv(np.reshape(cl_mat,(2,2,lmax-1)).transpose(2,0,1))
        dmat = []
        dmat.append([deriv_tt[:,2:],deriv_te[:,2:]])
        dmat.append([deriv_te[:,2:],deriv_ee[:,2:]])
        dmat = np.reshape(dmat,(2,2,len(cv_floats),lmax-1)).transpose(2,3,0,1)
        fish = []
        for i in range(len(cv_floats)):
            for j in range(len(cv_floats)):
                for k in range(len(ell[0,2:])):
                    fish.append(((2*ell[0][k]+1)*fsky/2)*np.trace(inv[k].dot(dmat[i][k].dot(inv[k]).dot(dmat[j][k]))))
        fish = np.reshape(fish,(len(cv_floats),len(cv_floats),len(ell[0][2:])))
        fish = np.sum(fish,axis = 2)
        return fish

    #planck -- return fisher(deriv,ell,cl_tt,40*np.pi/10800,7*np.pi/10800,0.65)
    #s4 -- return fisher(deriv,ell,cl_tt,1*np.pi/10800,3*np.pi/10800,0.60)
    return fisher(deriv,ell,cl_tt,0,3*np.pi/10800,0.50)

#params from table 1 of Planck 2015 results paper
#fisher_matrix = run_fisher(['output','lensing','omega_b','omega_cdm','H0','n_s','A_s','tau_reio'],['tCl,pCl,lCl','yes',0.02222,0.1199,67.26,0.9652,2.199e-9,0.078],[0.00023,0.0022,0.98,0.0062,0.016e-9,0.019])
fisher_matrix = run_fisher(['output','lensing','omega_b','omega_cdm','H0','n_s','A_s','tau_reio'],['tCl,pCl','no',0.0222,0.1197,67.31,0.9655,2.2e-9,0.06],[0.00022,0.0012,0.67,0.0097,0.022e-9,0.0006], run = False)
#print fisher_matrix


cov = np.linalg.inv(fisher_matrix)
cv = [0.02222,0.1199,67.26,0.9652,2.199e-9,0.078]
delta = [0.00023,0.0022,0.98,0.0062,0.016e-9,0.019]
wb = np.sqrt(cov[0][0])
wc = np.sqrt(cov[1][1])
h = np.sqrt(cov[2][2])
ns = np.sqrt(cov[3][3])
As = np.sqrt(cov[4][4])
tau = np.sqrt(cov[5][5])
print r'H0 : ',h
print r'ombh2 : ',wb
print r'omch2 : ',wc
print r'tau : ',tau
print r'As : ',As
print r'ns : ',ns

alpha_3 = 3.44
hold = [0,1,2,3,4,5]
comb = list(itertools.combinations(hold,2))
a2 = []
b2 = []
tan2th = []
for i in range(len(hold)):
    a2.append(((cov[comb[i][0]][comb[i][0]]+cov[comb[i][1]][comb[i][1]])/2.)+np.sqrt((np.square(cov[comb[i][0]][comb[i][0]]-cov[comb[i][1]][comb[i][1]])/4.)+np.square(cov[comb[i][0]][comb[i][1]])))
    b2.append(((cov[comb[i][0]][comb[i][0]]+cov[comb[i][1]][comb[i][1]])/2.)-np.sqrt((np.square(cov[comb[i][0]][comb[i][0]]-cov[comb[i][1]][comb[i][1]])/4.)+np.square(cov[comb[i][0]][comb[i][1]])))
    tan2th.append(2.*cov[comb[i][0]][comb[i][1]]/(np.square(cov[comb[i][0]][comb[i][0]])-np.square(cov[comb[i][1]][comb[i][1]])))
a = np.sqrt(a2)*alpha_3
b = np.sqrt(b2)*alpha_3
th = np.arctan(tan2th)/2.

def plot_ellipse(param1, param2):
    param_labels = ['\omega_b','\omega_{cdm}','H_0','n_s','A_s','\\tau']
    plt.figure()
    ax = plt.gca()
    ellipse = Ellipse(xy = (cv[param1],cv[param2]), width = a[param1], height = b[param2], angle = np.rad2deg(th[param1]), alpha = 0.4)
    ax.add_patch(ellipse)
    plt.xlim(cv[param1]-(alpha_3*delta[param1]),cv[param1]+(alpha_3*delta[param2]))
    plt.ylim(cv[param2]-(alpha_3*delta[param2]),cv[param2]+(alpha_3*delta[param2]))
    plt.xlabel(r'$%s$' %(param_labels[param1]))
    plt.ylabel(r'$%s$' %(param_labels[param2]))
    plt.savefig('plots/ellipse.pdf')
    subprocess.call('open plots/ellipse.pdf',shell = True)

#plot_ellipse(0,1)

