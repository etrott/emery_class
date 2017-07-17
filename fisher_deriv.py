from classy import Class
import numpy as np
import matplotlib.pyplot as plt
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

def run_fisher(param,mass,cv,delta,num,deg,run = 'False'):
    def run_CLASS(dict):
        cosmo = Class()
        cosmo.set(dict)
        cosmo.compute()
        output = cosmo.lensed_cl(2500)
        return output['ell'],output['tt'],output['te'],output['ee']

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
    print runs

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
        ell = np.reshape(ell,(len(cv_floats),num,2501))
        cl_tt = np.reshape(cl_tt,(len(cv_floats),num,2501))
        cl_te = np.reshape(cl_te,(len(cv_floats),num,2501))
        cl_ee = np.reshape(cl_ee,(len(cv_floats),num,2501))
        np.save('mat_ell',ell)
        np.save('mat_cl_tt',cl_tt)
        np.save('mat_cl_te',cl_te)
        np.save('mat_cl_ee',cl_ee)

    if run == 'True':
        save()

    def load():
        ell = np.load('mat_ell.npy')
        cl_tt = np.load('mat_cl_tt.npy')
        cl_te = np.load('mat_cl_te.npy')
        cl_ee = np.load('mat_cl_ee.npy')
        return ell,cl_tt,cl_te,cl_ee

    ell,cl_tt,cl_te,cl_ee = load()

    ell = np.transpose(ell,(0,2,1))
    cl_tt = np.transpose(cl_tt,(0,2,1))
    cl_te = np.transpose(cl_te,(0,2,1))
    cl_ee = np.transpose(cl_ee,(0,2,1))

    def poly_fit(x,y,deg):
        poly = []
        for i in range(len(y)):
            for j in range(len(y[i])):
                poly.append(np.polyfit(x[i],y[i][j],deg))
        poly = np.reshape(poly,(len(y),len(y[0]),deg+1))
        return poly

    def deriv(cl):
        poly_coeff = poly_fit(span,cl,deg)
#        print poly_coeff
        deriv_coeff = []
        for i in range(len(poly_coeff)):
            for j in range(len(poly_coeff[i])):
                deriv_coeff.append(np.poly1d(np.polyder(poly_coeff[i][j])))
#    deriv_coeff = np.reshape(deriv_coeff,(len(poly_coeff),len(poly_coeff[0]),len(poly_coeff[0][0])-1))
        deriv_coeff = np.reshape(deriv_coeff,(len(poly_coeff),len(poly_coeff[0]),1))
#        print deriv_coeff
        deriv = []
        for i in range(len(deriv_coeff)):
            for j in range(len(deriv_coeff[i])):
                deriv.append(deriv_coeff[i][j][0](span[i][12]))
        deriv = np.reshape(deriv,(len(deriv_coeff),len(deriv_coeff[0])))
        return deriv

    deriv = deriv(cl_tt)

    def plot(x,y):
        y_mod = x*(x+1.0)*y/(2*np.pi)
        param_labels = ['\omega_b','\omega_{cdm}','H_0']
        fig,ax = plt.subplots()
        for i in range(len(deriv)):
            ax.plot(x[0],y[i],label = r'$%s$' %(param_labels[i]))
#        plt.yscale('log')
        plt.xlabel(r'$l$')
        plt.ylabel(r'$\partial{C_l^{TT}}/\partial{p}$')
        plt.title(r'CMB Power Spectrum Derivatives')
        plt.legend(frameon = False, loc = 'upper right')
        plt.savefig('deriv.pdf')
        subprocess.call('open deriv.pdf',shell = True)

#    plot(ell[:,:,12],deriv)

    def fisher(deriv,ell,cl,npix,sigma,theta):
        sigma = theta/np.sqrt(8*np.log(2))
        ell = ell[:,:,12]
        cl = cl[:,:,12]
        npix = np.asarray(npix).reshape(1,1)
        sigma = np.asarray(sigma).reshape(1,1)
        theta = np.asarray(theta).reshape(1,1)
        matrix = []
#        n = 4*np.pi*np.square(sigma)/npix
#        w = np.exp(-np.square(ell)*np.square(theta)/(16*np.log(2)))
#        err = 2*(cl+n*np.power(w,-2))/(2*ell+1)
#        err = 2*(cl+n)/(2*ell+1)
        for i in range(len(deriv)):
            for j in range(len(deriv[i])):
                if np.isinf(err[i][j]) == True:
                    err[i][j] = 0
        for i in range(len(deriv)):
            for j in range(len(deriv)):
#                matrix.append(np.sum(deriv[i]*deriv[j]*err[i]))
                matrix.append(np.sum(deriv[i]*deriv[j]))
        matrix = np.reshape(matrix,(len(deriv),len(deriv)))
        return matrix

#    return fisher(deriv,ell,cl_tt,61,1,0.02217)


#fisher_matrix = run_fisher(['output','lensing','Omega_b','Omega_cdm','H0'],[1.0],['tCl,pCl,lCl','yes',0.02234,0.1189,67.8],[0.00023,0.0022,1.0],25,3,run = 'False')
fisher_matrix = run_fisher(['output','lensing','Omega_b','Omega_cdm','H0'],[1.0],['tCl,pCl,lCl','yes',0.0222,0.1197,67.31],[0.00022,0.0012,0.67],5,3,run = 'True')
print fisher_matrix
