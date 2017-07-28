bpars = {}

bpars['output'] = 'tCl pCl mPk'
bpars['gauge'] = 'synchronous'
bpars['l_max_scalars'] = 5000

bpars['h'] = 0.67556
bpars['T_cmb'] = 2.7255
bpars['Omega_b'] = 0.022032
bpars['N_ur'] = 3.046
bpars['Omega_k'] = 0.
bpars['z_reio'] = 11.357
bpars['A_s'] = 2.215e-9
bpars['n_s'] = 0.9619
bpars['alpha_s'] = 0.


bpars['reionization_exponent'] = 1.5
bpars['reionization_width'] = 0.5
bpars['helium_fullreio_redshift'] = 3.5
bpars['helium_fullreio_width'] = 0.5
bpars['lensing'] = 'no'


# bpars['l_max_tensors'] = 1000
# bpars['r'] = 0.
# ln10^{10}A_s = 3.0980
# tau_reio = 0.0925
# YHe = 'BBN'
# recombination = 'RECFAST'
# reio_parametrization = reio_camb

# 4) list of modes ('s' for scalars, 'v' for vectors, 't' for tensors). More than
#    one letter allowed, can be attached or separated by arbitrary characters;
#    letters can be small or capital.
#    (default: set to 's')
#modes = s,t

# 5) relevant only if you ask for 'tCl, lCl' and/or 'pCl, lCl': if you want the
#    spectrum of lensed Cls, enter a word containing the letter 'y' or 'Y'
#    (default: no lensed Cls)
# lensing = yes

# 6) which perturbations should be included in tensor calculations? write 'exact'
#    to include photons, ultra-relativistic species 'ur' and all non-cold dark
#    matter species 'ncdm'; write 'massless' to appriximate 'ncdm' as extra
#    relativistic species (good approximation if ncdm is still relativistic at
#    the time of recombination); write 'photons' to include only photons
#    (default: 'massless')
# tensor method =

# 7d) Do you want to write a table of background quantitites in a file? This will
#     include H, densities, Omegas, various cosmological distances, sound
#     horizon, etc., as a function of conformal time, proper time, scale factor.
#     File created if 'write background'  set to something containing the letter
#     'y' or 'Y', file written, otherwise not written (default: not written)
# write background = no

# 7e) Do you want to write a table of thermodynamics quantitites in a file? File
#     is created if 'write thermodynamics' is set to something containing the
#     letter 'y' or 'Y'. (default: not written)
# write thermodynamics = no

# 7f) Do you want to write a table of perturbations to files for certain
#     wavenumbers k? Dimension of k is 1/Mpc. The actual wave numbers are chosen
#     such that they are as close as possible to the requested k-values.
# k_output_values = #0.01, 0.1, 0.0001
