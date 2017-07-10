from classy import Class
import numpy as np

ell = np.arange(2501).reshape(1,2501)

def run_fisher(param_list,masses,central_values,fin_epsilon_definite,run = 'False'):
    def run_CLASS(dict):
        cosmo = Class()
        cosmo.set(dict)
        cosmo.compute()
        output = cosmo.lensed_cl(2500)
        return output['tt']
    num_cv = []
    num_ep = []
    fin_epsilon = []
    for i in range(len(param_list)):
        if isinstance(central_values[i],basestring):
            fin_epsilon.append(fin_epsilon_definite[i])
        else:
            num_cv.append(central_values[i])
            num_ep.append(fin_epsilon_definite[i])
            fin_epsilon.append(fin_epsilon_definite[i]/central_values[i])
    num_cv = np.asarray(num_cv).reshape(len(num_cv),1)
    num_ep = np.asarray(num_ep).reshape(len(num_cv),1)
    param_values = []
    for i in range(len(param_list)):
        if isinstance(central_values[i],basestring):
            param_values.append((central_values[i],central_values[i]))
        else:
            param_values.append((central_values[i]*(1.0-fin_epsilon[i]),central_values[i]*(1.0+fin_epsilon[i])))
    dict_low = {param_list[i]:param_values[i][0] for i in range(len(param_list))}
    dict_high = {param_list[i]:param_values[i][1] for i in range(len(param_list))}
    def save(dict_low,dict_high):
        cl_low = []
        cl_high = []
        cl_low.append(run_CLASS(dict_low))
        cl_high.append(run_CLASS(dict_high))
        cl_low = np.reshape(cl_low,(len(num_cv),2501))
        cl_high = np.reshape(cl_high,(len(num_cv),2501))
        np.save('cl_low',cl_low)
        np.save('cl_high',cl_high)

    if run == 'True':
        save(dict_low,dict_high)

    def load():
        cl_low = np.load('cl_low.npy')
        cl_high = np.load('cl_high.npy')
        cl_low = np.reshape(cl_low,(len(num_cv),len(masses),len(ell[0])))
        cl_high = np.reshape(cl_high,(len(num_cv),len(masses),len(ell[0])))
        return cl_low,cl_high
    cl_low,cl_high = load()

    def fin_deriv(cl_high,cl_low,central_value,step_size):
        param_deriv = (cl_high-cl_low)/(central_value*step_size)
        return param_deriv
    deriv = fin_deriv(cl_high,cl_low,num_cv,num_ep)
    def fisher(deriv):
        matrix = []
        for i in range(len(deriv)):
            for j in range(len(deriv)):
                matrix.append(np.sum(deriv[i]*deriv[j]))
        matrix = np.reshape(matrix,(len(num_cv),len(num_cv)))
        return matrix
    return fisher(deriv)

fisher_matrix = run_fisher(['output','lensing','Omega_b','Omega_cdm','H0'],[1.0],['tCl,pCl,lCl','yes',0.022032,0.12038,67.556],['string','string',0.00125,0.005,0.08],run = 'False')
print fisher_matrix
