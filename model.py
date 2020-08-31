#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Handy MCMC scripts.

Test for the different fit method (mcmc, ptmcmc, minimizer).

Author:
    Dimitri Misiak (misiak@ipnl.in2p3.fr)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

# importing ethem package
from import_package import custom_import
custom_import()
import ethem as eth

from config_ethem import syst, get_eval_dict
from raw_data import get_raw_data

def initiate_model():
    """ Select the parameters for the model.

    Return
    ======
    param_top: tuple
        Symbols of the primary parameters of the models.

    param_bot: tuple
        Symbols of the secondary parameters of the models.

    sys_noise_full: fun
        Simulation function from ethem.
    """
    cryo = syst.Thermostat_b
    capa = syst.Capacitor_f
    elntd = syst.Resistor_ntd
    leak = syst.ThermalLink_leak
    epcoup = syst.ThermalLink_ep
#    glue = syst.ThermalLink_glue
    
    # reference bath
#    ref_bath = syst.Capacitor_f
    
    # "primary" parameters
    param_top = (cryo.temperature,
                 capa.current,
                 )

    # "secondary" parameters
    param_bot = (elntd.R0,
                 elntd.T0,
                 leak.cond_alpha,
                 epcoup.cond_alpha,
                 )

    param_full = param_top + param_bot

    # evaluation dictionnary
    evad = get_eval_dict()

    ss_point = eth.solve_sse_param(syst, param_full, evad)

    return param_top, param_bot, ss_point

### dirty !!! should be better with classes
t_array, iv_list = get_raw_data()

##########################
#WORK IN PROGRESS
##########################
def make_model_data(param):
    """ Return the model array from the given set of parameters.
    """
    param_top, param_bot, ss_point = initiate_model()

    model_iv_list = list()
    
    for i,temp in enumerate(t_array):

        i_array, data_iv = iv_list[i]

        model_iv = list()
        for curr in i_array:
            
            p_top = (temp, curr)
            p_full = p_top + tuple(param)
            
            model_iv.append(ss_point(p_full).x[-1])
            
        model_iv_array = np.vstack((i_array, model_iv))
        model_iv_list.append(model_iv_array)

    return model_iv_list


def make_fake_data(model_data, sigma_coeff=1.):
    """ Add a gaussian noise depending on the model_data.
    """
    fake_data = list()
    for md in model_data:

        freq, y_data = md

        sigma = y_data * sigma_coeff

        blob = np.random.normal(0,sigma,sigma.shape)

        noisy_data = y_data + blob

        fake_data.append(np.vstack((freq, noisy_data)))

    return fake_data


if __name__ == "__main__":

    param_top, param_bot, ss_point = initiate_model()
    evad = get_eval_dict()
    p0 = [evad[p] for p in param_bot]

    model0 = make_model_data(p0)
    
    fake0 = make_fake_data(model0, sigma_coeff=0.1)

    t_array, iv_list = get_raw_data()

    cmap = plt.get_cmap('magma')

    plt.close('all')
    plt.figure('noise dict plot', figsize=(10,7))

    for i,temp in enumerate(t_array):
        
#        i_array, data_iv = iv_list[i]
        i_array, fake_iv = fake0[i]
        _, model_iv = model0[i]
        
        color = cmap((i+1.)/(len(t_array)+1))

        plt.loglog(
                i_array, fake_iv,
                label='{} K'.format(temp),
                color=color,
                lw=0.5,
                marker='.',
        )

        plt.loglog(
                i_array, model_iv,
                label='{} K'.format(temp),
                color=color,
                lw=1.5,
                path_effects=[pe.Stroke(linewidth=2, foreground='k'),
                              pe.Normal()],
        )

    plt.legend()
    plt.grid(True)
