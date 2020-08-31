#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Plotting the experimental data
@author: misiak
"""

import matplotlib.pyplot as plt
import numpy as np

from import_package import custom_import


def get_raw_data():
    """ Return the raw noise data. For a different set of data, edit this
    function.

    Return
    ======
    name_array: array of stream names.

    temp_array: array of temperature.

    res_array: array of NTD resistance.

    xy_data_list: list of freq,lpsd arrays.
    """

    data_dir, output_dir = custom_import()

    vi_path = ''.join((data_dir, '/vi_run53_red50.csv'))

    i_data, v_data, t_data = np.loadtxt(vi_path, comments='#', unpack=True)

    t_array = np.unique(t_data)

    iv_list = list()
    for t in t_array:
        index_array = np.where(t_data == t)[0]
        i_array = i_data[index_array]
        v_array = v_data[index_array]
      
        iv_list.append( np.vstack((i_array, v_array)) )

    return t_array, iv_list

if __name__ == '__main__':
    
    t_array, iv_list = get_raw_data()
    
    plt.figure(figsize=(8,5))
    
    for t, iv_array in zip(t_array, iv_list):
        i_array, v_array = iv_array
        plt.errorbar(i_array, v_array, yerr=v_array*0.1,
                     lw=1., ls='-',
                     label='{0:.1f} mK'.format(t*1e3))
    
    plt.grid()
    plt.xlabel('Voltage [V]')
    plt.ylabel('Current [A]')
    plt.ylim(1e-8, 1e-2)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout(rect=(0., 0., 0.8, 1.))
    
    
    print('Done')