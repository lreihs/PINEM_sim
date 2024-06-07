import numpy as np
import scipy.constants as const

params = []

f = open('../docs/input_params', 'r')
for line in f:
    save_entry = False
    for entry in line.split():
        if save_entry == True:
            params.append(float(entry))
        if list(entry)[-1] == '=':
            save_entry = True

[logN_samp,
 z_max_mm,
 t_max_stdev,
 E_el_kev, 
 E_el_fwhm_ev, 
 pulse_duration_fs,
 lambda_nm, 
 field_nm, 
 min_field_nm, 
 max_field_nm, 
 N_spec, 
 N_steps, 
 dt] = params

N_samp = int(10**logN_samp)
N_spec = int(N_spec)

hbar = const.hbar
e = const.elementary_charge
pi = const.pi
eV = const.electron_volt
c = const.speed_of_light

z_max = z_max_mm*10**-3
E_el = E_el_kev*10**3*eV
E_el_fwhm = E_el_fwhm_ev*eV
pulse_duration = pulse_duration_fs*10**-15
field = field_nm*10**9
min_field = min_field_nm*10**9
max_field = max_field_nm*10**9
lambda_opt = lambda_nm*10**-9
w_opt = 2*pi*c/lambda_opt 

z = np.linspace(-z_max/2, z_max/2, N_samp)
t = np.linspace(-pulse_duration*t_max_stdev/2, pulse_duration*t_max_stdev/2, N_samp)
w = 2*pi*np.fft.fftshift(np.fft.fftfreq(N_samp, d = pulse_duration*t_max_stdev/N_samp))
photon_order = w/w_opt 
g = 49*10**-9*field
