from classy import Class
import numpy as np
from numpy import linalg

master = {'run': False,
'dmeff_mass': 0.001,
'fname': '1MeV',
'survey': 'planck', 
'polarization': True,
'param':['output','lensing','omega_b','omega_dmeff','H0','n_s','A_s','tau_reio','cc_dmeff_p'],
'cv':['tCl,pCl','no',0.0222,0.1197,67.31,0.9655,2.2e-9,0.06,0.],
'delta':[0.00022,0.0012,0.67,0.0097,0.022e-9,0.0006,5e8]}

def run_CLASS(dict):
    cosmo = Class()
    cosmo.set(dict)
    cosmo.compute()
    output = cosmo.raw_cl(lmax)
    #output = cosmo.lensed_cl(lmax)
    return output['ell'],output['tt'],output['te'],output['ee']

def sep_dict(param,cv):
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
    return cv_strings,param_strings,cv_floats,param_floats,mid

def construct_run(cv_floats,delta,num):
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
                span[i] = np.linspace(cv_floats[i],high[i][i],num)
    runs = np.reshape([cv_floats]*len(cv_floats)*num,(len(cv_floats),num,len(cv_floats)))
    for i in range(len(runs)):
        runs[i][:,i] = span[i]
    return runs,span

def save_pcc(fname,dmeff_mass,runs,cv_strings,param_strings,cv_floats,param_floats):
    ell = []
    cl_tt = []
    cl_te = []
    cl_ee = []
    for i in range(len(runs)):
        if i == 6:
            for j in range(len(runs[i])):
                dict = {param_floats[k]:runs[i][j][k] for k in range(len(param_floats))}
                string_dict = {param_strings[k]:cv_strings[k] for k in range(len(param_strings))}
                dmeff_dict = {'m_dmeff':dmeff_mass,'cc_dmeff_op':1, 'cc_dmeff_num':1, 'cc_dmeff_n':0, 'cc_dmeff_qm2':0, 'cc_dmeff_scale': 1, 'omega_cdm': 0.0, 'spin_dmeff':0., 'use_helium_dmeff':'yes', 'use_temperature_dmeff':'yes'}
                dict.update(string_dict)
                dict.update(dmeff_dict)
                ell.append(run_CLASS(dict)[0])
                cl_tt.append(run_CLASS(dict)[1])
                cl_te.append(run_CLASS(dict)[2])
                cl_ee.append(run_CLASS(dict)[3])
    ell = np.reshape(ell,(1,num,lmax+1))
    cl_tt = np.reshape(cl_tt,(1,num,lmax+1))
    cl_te = np.reshape(cl_te,(1,num,lmax+1))
    cl_ee = np.reshape(cl_ee,(1,num,lmax+1))
    np.save('data/ell_%s'%(fname),ell)
    np.save('data/cl_tt_%s'%(fname),cl_tt)
    np.save('data/cl_te_%s'%(fname),cl_te)
    np.save('data/cl_ee_%s'%(fname),cl_ee)

def save(fname,dmeff_mass,runs,cv_strings,param_strings,cv_floats,param_floats):
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
            dmeff_dict = {'m_dmeff':dmeff_mass,'cc_dmeff_op':1, 'cc_dmeff_num':1, 'cc_dmeff_n':0, 'cc_dmeff_qm2':0, 'cc_dmeff_scale': 1, 'omega_cdm': 0.0, 'spin_dmeff':0., 'use_helium_dmeff':'yes', 'use_temperature_dmeff':'yes'}
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

def load(fname):
    ell = np.load('data/ell_%s.npy'%(fname))
    cl_tt = np.load('data/cl_tt_%s.npy'%(fname))
    cl_te = np.load('data/cl_te_%s.npy'%(fname))
    cl_ee = np.load('data/cl_ee_%s.npy'%(fname))
    return ell,cl_tt,cl_te,cl_ee

def temp(ell,cl_tt,cl_te,cl_ee):
    scale = np.asarray((2.7255**2)*(10**12)).reshape(1,1,1)
    ell = np.transpose(ell,(0,2,1))
    cl_tt = np.transpose(cl_tt,(0,2,1))
    cl_te = np.transpose(cl_te,(0,2,1))
    cl_ee = np.transpose(cl_ee,(0,2,1))
    cl_tt = cl_tt*scale
    cl_te = cl_te*scale
    cl_ee = cl_ee*scale
    return ell,cl_tt,cl_te,cl_ee

def poly_fit(x,y,deg):
    poly = []
    for i in range(len(y)):
        for j in range(len(y[i])):
            poly.append(np.polyfit(x[i],y[i][j],deg))
    poly = np.reshape(poly,(len(y),len(y[0]),deg+1))
    return poly

def deriv(span,cl,cv_floats):
    poly_coeff = poly_fit(span,cl,deg)
    deriv_coeff = []
    for i in range(len(poly_coeff)):
        for j in range(len(poly_coeff[i])):
            deriv_coeff.append(np.poly1d(np.polyder(poly_coeff[i][j])))
    deriv_coeff = np.reshape(deriv_coeff,(len(poly_coeff),len(poly_coeff[0]),1))
    deriv = []
    for i in range(len(deriv_coeff)):
        for j in range(len(deriv_coeff[i])):
            deriv.append(deriv_coeff[i][j][0](cv_floats[i]))
    deriv = np.reshape(deriv,(len(deriv_coeff),len(deriv_coeff[0])))
    return deriv

def fisher(ell,cl_tt,cl_te,cl_ee,deriv_tt,deriv_te,deriv_ee,cv_floats,param_floats,s,theta,fsky,polarization = True):
    ell = ell[0,:,mid][2:]
    a = cl_tt[0,:,mid][2:]
    b = cl_te[0,:,mid][2:]
    c = cl_ee[0,:,mid][2:]
    s = np.asarray(s)
    theta = np.asarray(theta)
    n = np.square(s)*np.exp(ell*(ell+1.0)*np.square(theta)/(8.0*np.log(2.0)))
    cl_mat = []
    dmat = []
    if polarization == True:
        cl_mat.append([a+n,b])
        cl_mat.append([b,c+(2*n)])
        cl_mat = np.transpose(cl_mat,(2,0,1))
        inv = np.linalg.inv(cl_mat)
        dmat.append([deriv_tt[:,2:],deriv_te[:,2:]])
        dmat.append([deriv_te[:,2:],deriv_ee[:,2:]])    
        dmat = np.transpose(dmat,(2,3,0,1))
    if polarization == False:
        cl_mat.append([a+n])
        cl_mat = np.transpose(cl_mat,(2,0,1))
        inv = np.linalg.inv(cl_mat)
        dmat.append([deriv_tt[:,2:]])
        dmat = np.transpose(dmat,(2,3,0,1))
    fish = []
    for i in range(len(cv_floats)):
        for j in range(len(cv_floats)):
            for k in range(len(ell)):
                fish.append(((2.*ell[k]+1.)*fsky/2.)*np.trace(inv[k].dot(dmat[i][k].dot(inv[k]).dot(dmat[j][k]))))
    fish = np.reshape(fish,(len(cv_floats),len(cv_floats),len(ell)))
    fish = np.sum(fish,axis = 2)
    return fish

def noise(ell,s,theta):
    ell = ell[0,:,mid][2:]
    s = np.asarray(s)
    theta = np.asarray(theta)
    n = np.square(s)*np.exp(ell*(ell+1.0)*np.square(theta)/(8.0*np.log(2.0)))
    return n

def output_constraints(fisher_matrix,param):
    cov = np.linalg.inv(fisher_matrix)
    for i in range(len(param)):
        constraint = np.sqrt(cov[i][i])
        print '%s : '%(param[i]), constraint

def cross_section(fisher_matrix,dmeff_mass,param_floats,proton_mass = 0.938272,v = 246.):
    cov = np.linalg.inv(fisher_matrix)
    mu = (dmeff_mass*proton_mass)/(dmeff_mass+proton_mass)
    for i in range(len(param_floats)):
        if param_floats[i] == 'cc_dmeff_p':
            sigma = ((cov[i][i]*np.square(mu))/(np.pi*np.power(v,4.)))/np.square(2.5*10**13)
            print 'p_cc : ',sigma
"""
if master['polarization'] == True:
    print master['survey'],master['fname'],'with polarization'
if master['polarization'] == False:
    print master['survey'],master['fname'],'w/o polarization'
"""
#param = ['output','lensing','omega_b','omega_dmeff','H0','n_s','A_s','tau_reio','cc_dmeff_p']
#cv = ['tCl,pCl','no',0.0222,0.1197,67.31,0.9655,2.2e-9,0.06,0.]
#delta = [0.00022,0.0012,0.67,0.0097,0.022e-9,0.0006,5e8]
num = 25
deg = 3
lmax = 2500
mid = (num-1)/2

def forecast(master):
    cv_strings,param_strings,cv_floats,param_floats,mid = sep_dict(master['param'],master['cv'])
    runs,span = construct_run(cv_floats,master['delta'],num)
    if master['run'] == True:
        save(master['fname'],master['dmeff_mass'],runs,cv_strings,param_strings,cv_floats,param_floats)
    ell,cl_tt,cl_te,cl_ee = load(master['fname'])
    ell,cl_tt,cl_te,cl_ee = temp(ell,cl_tt,cl_te,cl_ee)
    deriv_tt,deriv_te,deriv_ee = deriv(span,cl_tt,cv_floats),deriv(span,cl_te,cv_floats),deriv(span,cl_ee,cv_floats)
    if master['survey'] == 'planck':
        fisher_matrix = fisher(ell,cl_tt,cl_te,cl_ee,deriv_tt,deriv_te,deriv_ee,cv_floats,param_floats,40.*np.pi/10800.,7.*np.pi/10800.,0.5,polarization = master['polarization'])
    if master['survey'] == 's4':
        fisher_matrix = fisher(ell,cl_tt,cl_te,cl_ee,deriv_tt,deriv_te,deriv_ee,cv_floats,param_floats,1.*np.pi/10800.,3.*np.pi/10800.,0.6,polarization = master['polarization'])
    if master['survey'] == 'cv_limited':
        fisher_matrix = fisher(ell,cl_tt,cl_te,cl_ee,deriv_tt,deriv_te,deriv_ee,cv_floats,param_floats,0.,3.*np.pi/10800.,0.5,polarization = master['polarization'])
    output_constraints(fisher_matrix,['omega_b','omega_dmeff','H0','n_s','A_s','tau_reio','cc_dmeff_p'])
    cross_section(fisher_matrix,master['dmeff_mass'],param_floats)

def output_cl(master):
    cv_strings,param_strings,cv_floats,param_floats,mid = sep_dict(master['param'],master['cv'])
    runs,span = construct_run(cv_floats,master['delta'],num)
    #if master['run'] == True:
    #    save(master['fname'],master['dmeff_mass'],runs,cv_strings,param_strings,cv_floats,param_floats)
    ell,cl_tt,cl_te,cl_ee = load(master['fname'])
    ell,cl_tt,cl_te,cl_ee = temp(ell,cl_tt,cl_te,cl_ee)
    return ell,cl_tt,cl_te,cl_ee

def output_deriv(master):
    cv_strings,param_strings,cv_floats,param_floats,mid = sep_dict(master['param'],master['cv'])
    runs,span = construct_run(cv_floats,master['delta'],num)
    #if master['run'] == True:
    #    save(master['fname'],master['dmeff_mass'],runs,cv_strings,param_strings,cv_floats,param_floats)
    ell,cl_tt,cl_te,cl_ee = load(master['fname'])
    ell,cl_tt,cl_te,cl_ee = temp(ell,cl_tt,cl_te,cl_ee)
    deriv_tt,deriv_te,deriv_ee = deriv(span,cl_tt,cv_floats),deriv(span,cl_te,cv_floats),deriv(span,cl_ee,cv_floats)
    return deriv_tt,deriv_te,deriv_ee


#if __name__ == "__main__":