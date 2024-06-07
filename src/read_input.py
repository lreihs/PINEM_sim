import numpy as np
import scipy.constants as const

params = []

f = open('../docs/input_params', 'r') # change directory name if you're not Lennard
for line in f:
    save_entry = False
    for entry in line.split():
        if save_entry == True:
            params.append(float(entry))
        if list(entry)[-1] == ':':
            save_entry = True

[atomic_coordinates, 
 N_samp,
 z_max_mm,
 t_max_stdev,
 E_el_kev, 
 E_el_fwhm_ev, 
 pulse_duration_fs, 
 first_order_chirp_coeff_fs, 
 lambda_nm, 
 field_nm, 
 min_field_nm, 
 max_field_nm, 
 N_spec, 
 N_steps, 
 dt] = params
