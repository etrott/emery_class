import os,os.path,shutil,sys
import numpy as np
import pickle
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except:
    print 'Warning: did not import matplotlib, plotting functinos will not work.'

from classy import Class
from classy import CosmoSevereError
from baseline_inputs import bpars

BPARS = bpars.copy()
OMEGA0_DM = 0.12038
PATH = '/Users/verag/Research/Repositories/DMeff-CMB/class_tests/'

###############################################################################
def compare_cl_dicts(cl1, cl2, atol=1e-20,rtol=1e-3):
    """
    compares two cl dictionaries in the format as classy output.
    """
    
    for k,v in cl1.iteritems():
        if k in cl2.keys():
            if np.allclose(v,cl2[k],atol=atol,rtol=rtol):
                pass
            else:
                return False
        else:
            return False
    return True

###############################################################################
def make_baseline(model,
                  baseline_name='cl_baseline',
                  m_dmeff=1.,
                  spin_dmeff=0.5,
                  cc_dmeff_qm2=0,
                  operators=1,
                  cc_dmeff_p=1000,
                  cc_dmeff_n=500):
    """
    computes and saves baseline files

    Input:

    model = either 'cdm' or 'dmeff'
    m_dmeff = mass of DM particle in GeV (for model=dmeff)
    spin_dmeff = DM spine (for model=dmeff)
    cc_dmeff_qm2 = additional momentum dependence of the interaction (usually 0; for model=dmeff)
    operators = must be at least one integer, or a list of up to 15 integers <=15.
    cc_dmeff_p = coupling const with protons; one float, or a list of up to 15 floats, must match operators size.
    cc_dmeff_n = coupling const with neutrons; one float, or a list of up to 15 floats, must match operators size.
    """
    operators = np.atleast_1d(operators)
    cc_dmeff_p = np.atleast_1d(cc_dmeff_p)
    cc_dmeff_n = np.atleast_1d(cc_dmeff_n)
    if model=='dmeff':
        if (len(operators) != len(cc_dmeff_p)) or (len(operators) != len(cc_dmeff_n)):
            raise ValueError('lengths of operators, cc_dmeff_p, and cc_dmeff_n do not match!')

    for i,op in enumerate(operators):
        print('\n Calculating Cls for {}'.format(model))

        basefile_name = '{}_{}'.format(baseline_name, model)
        if model=='dmeff':
            basefile_name += '_op{}_{}GeV_spin{}_qm2{}_ccp{}_ccn{}'.format(op,m_dmeff,spin_dmeff,cc_dmeff_qm2,cc_dmeff_p[i],cc_dmeff_n[i])
            print('operator {}...\n').format(op)
        basefile_name += '.pkl'

        kwargs = {}
        if model=='dmeff':
            kwargs['m_dmeff'] = m_dmeff
            kwargs['spin_dmeff'] = spin_dmeff
            kwargs['cc_dmeff_op'] = op
            kwargs['cc_dmeff_qm2'] = cc_dmeff_qm2
            kwargs['cc_dmeff_p'] = cc_dmeff_p[i]
            kwargs['cc_dmeff_n'] = cc_dmeff_n[i]

        cl = compute_cls(model=model,**kwargs)
        pickle.dump( cl, open( basefile_name, "wb" ) )

    
###############################################################################
def compute_cls(model='cdm',**kwargs):
    """
    model can be cdm or dmeff.
    """

    par = BPARS.copy()
    if model=='dmeff':
        par['omega_cdm'] = 0.
        par['omega_dmeff'] = OMEGA0_DM
        par['cc_dmeff_scale'] = 1
        par['cc_dmeff_num'] = 1
        par['use_helium_dmeff']      = 'yes'
        par['use_temperature_dmeff'] = 'no'
    elif model=='cdm':
        par['omega_dmeff'] = 0.
        par['omega_cdm'] = OMEGA0_DM
        
    for k in kwargs.keys():
        par[k] = kwargs[k]
    
    cosmo = Class() # Create CLASS instance
    cosmo.set(par)  # Set CLASS parameters
    cosmo.compute() # CLASS computation

    cl = cosmo.raw_cl()
    
    # Free memory (equivalent of struct_free() in `main` of CLASS).
    # This is essential when looping over different cosmologies.
    cosmo.struct_cleanup()

    # Clean up parameters, though it is not nedded if looping over
    # different values of the same set of parameters.
    cosmo.empty()

    return cl

###############################################################################
def plot_cls(Cls,plot_ratio=True,
             labels=None, plotname='plot.pdf',
             styles=None,saveplot=True):
    """
    Cls is a list of cl objects (dicts) that classy returns
    """
    if styles is None:
        styles = ['-']*len(Cls)
    if labels is None:
        labels = ['']*len(Cls)

    plt.figure()
    
    plotted = False
    if plot_ratio:
        plt.ylabel(r'$\delta C_\ell/C_\ell$')
        cl1 = Cls[0]
        cl2 = Cls[1]
        l = cl1['ell']
        for k in cl1.keys():
            if k!='ell':
                if ( (cl1[k][2:]==0.).sum() > 0) or ( (cl2[k][2:]==0.).sum() > 0):
                    print 'for {} there are some cls==0, so skipping plotting...'.format(k)
                else:
                    if ( (( cl1[k][2:] - cl2[k][2:] )==0.).sum() > 0. ):
                        print 'for {} the cls are the same, so skipping plotting...'.format(k)
                    else:
                        plt.loglog(l[2:], np.abs( ( cl1[k][2:] - cl2[k][2:] ) / cl1[k][2:] ),label=k)
                        plotted = True
                   
    else:
        plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$')
        for i,cl in enumerate(Cls):
            l = cl['ell']
            lf = l * (l + 1.) / 2./np.pi
            for k,v in cl.iteritems():
                if k!='ell':
                    plt.loglog(l, v*lf,styles[i],label=k + labels[i])
        
   

    plt.xlabel(r'$\ell$')
    
    if plotted:
        plt.legend(loc='lower left')
    
    if saveplot:
        plt.savefig(plotname)
    plt.close()
            
###############################################################################
def check_cl_match(baseline_file, model,
                   makeplot=True, saveplot=True,
                   atol=1e-20,rtol=1e-3,
                   **kwargs):
    """
    checks that current version of class
    returns the same cls as those pre-computed and stored in baseline_file
    """
    cl = compute_cls(model=model,**kwargs)
    l = cl['ell']
    
    clb = pickle.load( open( baseline_file, "rb" ) )
    lb = clb['ell']

    
    if makeplot:
        try:
            plot_cls([clb,cl],labels=[' (baseline)', ' (new)'],styles=['-','--'],
                     plotname=baseline_file[:-4] + '__test.pdf',
                     saveplot=saveplot)
        except:
            #raise
            print('(class_tests: failed to plot)\n')

    print cl['tt'],clb['tt']
    assert compare_cl_dicts(cl, clb, atol=atol, rtol=rtol)
    
###############################################################################
def test_cdm(baseline_file='cl_baseline_cdm.pkl',atol=1e-20,rtol=1e-3):
    """
    main testing function to match cl baseline with cdm model.
    """
    print('\nTesting Cls for cdm...\n')
    check_cl_match(baseline_file,'cdm',atol=atol,rtol=rtol)


###############################################################################
def test_dmeff(m_dmeff=1.,
                spin_dmeff=0.5,
                cc_dmeff_qm2=0,
                atol=1e-20,rtol=1e-3):
    """
    main testing function to match cl baseline with dmeff model, for all 15 operators.

    Input:

    m_dmeff = mass of DM particle in GeV (for model=dmeff)
    spin_dmeff = DM spine (for model=dmeff)
    cc_dmeff_op = interaction operators to turn on (for model=dmeff)
    cc_dmeff_qm2 = additional momentum dependence of the interaction (usually 0; for model=dmeff)
    """

    operators = [1]
    cc_dmeff_p = [1000]
    cc_dmeff_n = [500]
    
    for i,op in enumerate(operators):
        print('\nTesting Cls for dmeff, operator {}...\n'.format(op))

        baseline_file = 'cl_baseline_dmeff_op{}_{}GeV_spin{}_qm2{}_ccp{}_ccn{}.pkl'.format(op,m_dmeff,spin_dmeff,cc_dmeff_qm2,cc_dmeff_p[i],cc_dmeff_n[i])
        
        kwargs = {}
        kwargs = BPARS.copy()
        kwargs['m_dmeff'] = m_dmeff
        kwargs['spin_dmeff'] = spin_dmeff
        kwargs['cc_dmeff_op'] = op
        kwargs['cc_dmeff_qm2'] = cc_dmeff_qm2
        kwargs['cc_dmeff_p'] = cc_dmeff_p[i]
        kwargs['cc_dmeff_n'] = cc_dmeff_n[i]

        check_cl_match(baseline_file, 'dmeff',atol=atol,rtol=rtol, **kwargs)

###############################################################################
def test_cc0_limit(m_dmeff=1.,
                   spin_dmeff=0.5,
                   cc_dmeff_op=1,
                   cc_dmeff_qm2=0,
                   cc_dmeff_p=0.,
                   cc_dmeff_n=0.,
                   atol=1e-20,
                   rtol=1e-3,
                   makeplot=False):
    """
    checks that cc->0 recovers cdm result
    """
    
    kwargs = {}
    kwargs['m_dmeff'] = m_dmeff
    kwargs['spin_dmeff'] = spin_dmeff
    kwargs['cc_dmeff_op'] = cc_dmeff_op
    kwargs['cc_dmeff_qm2'] = cc_dmeff_qm2
    kwargs['cc_dmeff_p'] = cc_dmeff_p
    kwargs['cc_dmeff_n'] = cc_dmeff_n

    cl_cdm = compute_cls(model='cdm')
    cl_cc0 = compute_cls(model='dmeff', **kwargs)
    if makeplot:
        plot_cls([cl_cdm,cl_cc0],labels=[' (cdm)', ' (cc0)'],styles=['-','--'],
                     plotname='test_cc0_plot.pdf',
                     saveplot=True)
    #return cl_cdm,cl_cc0
    assert compare_cl_dicts(cl_cdm, cl_cc0, atol=atol,rtol=rtol)

######################################################
######################################################
if __name__=="__main__":
    test_cdm()
    #test_dmeff()
    
# kwargs = {}
# kwargs['m_dmeff'] = 1
# kwargs['spin_dmeff'] = 0.5
# kwargs['cc_dmeff_op'] = 1
# kwargs['cc_dmeff_qm2'] = 0
# kwargs['cc_dmeff_p'] = 0
# kwargs['cc_dmeff_n'] = 0
# kwargs['omega_cdm'] = 0.
# kwargs['omega_dmeff'] = 0.1
# kwargs['cc_dmeff_scale'] = 1
# kwargs['cc_dmeff_num'] = 1
# kwargs['use_helium_dmeff']      = 'yes'
# kwargs['use_temperature_dmeff'] = 'no'

# kwargs={'cc_dmeff_n': 0,
#  'cc_dmeff_num': 1,
#  'cc_dmeff_op': 1,
#  'cc_dmeff_p': 0,
#  'cc_dmeff_qm2': 0,
#  'cc_dmeff_scale': 1,
#  'm_dmeff': 1,
#  'omega_cdm': 0.0,
#  'omega_dmeff': 0.12038,
#  'spin_dmeff': 0.5,
#  'use_helium_dmeff': 'yes',
#  'use_temperature_dmeff': 'no'}
