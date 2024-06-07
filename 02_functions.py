import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import sys
from scipy.optimize import fsolve

### 2.1 Importing parameters from input file ###

params = []

f = open('/Users/lennardreihs/EPFL/Master/MA2/LUMES/Code/Current/01_params', 'r') # change directory name if you're not Lennard
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

N_samp = int(N_samp)
N_spec = int(N_spec)

if atomic_coordinates == True:
    hbar = 1
    e = 1
else:
    hbar = const.hbar
    e = const.elementary_charge

### 2.2 Additional paramters ###

def k_(E):
    return hbar**-1 * np.sqrt((2*const.electron_mass*E)*(1+E/(2*const.electron_mass*const.speed_of_light**2)))

def E_array(k):
    # k must be array
    return fsolve(lambda E0: E0**2 + 2*E0*const.electron_mass*const.speed_of_light**2 - hbar**2*k**2*const.speed_of_light**2, np.full(len(k), const.electron_volt))

def E_single(k):
    # k must be float
    return fsolve(lambda E0: E0**2 + 2*E0*const.electron_mass*const.speed_of_light**2 - hbar**2*k**2*const.speed_of_light**2, const.electron_volt)[0]

def dk_dE(E):
    return hbar**-1 * (const.electron_mass + E/const.speed_of_light**2) * (2*const.electron_mass*E + E**2/const.speed_of_light**2)**-0.5

def dE_dk(k):
    return hbar**2*k*(const.electron_mass**2 + hbar**2*k**2/const.speed_of_light**2)**-0.5

z_max = z_max_mm*10**-3
E_el = E_el_kev*const.electron_volt*10**3 # electron energy E0 [J]
E_el_fwhm = E_el_fwhm_ev*const.electron_volt # FWHM of inital electron energy spectrum [J]
k_el = k_(E_el) # frequency of plane wave component of unperturbed electron state [1/m]
k_el_fwhm = dk_dE(E_el_fwhm)*E_el_fwhm
v_el = (hbar * k_el / const.electron_mass) / (1 + E_el / (const.electron_mass*const.speed_of_light**2)) # electron z-velocity [m/s]
# lorentz_factor = 1/np.sqrt(1-v_el**2/const.speed_of_light**2) # Lorentz factor for electron with velocity v_el (defined above)
pulse_duration = pulse_duration_fs*10**(-15) # electron pulse duration (FWHM) [s]
first_order_chirp_coeff = first_order_chirp_coeff_fs * 10**-15 / const.electron_volt # first order chirp coefficnet of electron pulse [s/J]
field = field_nm*10**9 # field strength of optical field [V/m]
min_field = min_field_nm*10**9 # minimum field strength of optical field for field dependent spectra [V/m]
max_field = max_field_nm*10**9 # maximum field strength of optical field for field dependent spectra [V/m]
w_opt = 2 * const.pi * const.speed_of_light / (lambda_nm*10**-9) # angular frequency of light field [Hz]

z = np.linspace(-z_max/2, z_max/2, N_samp)
k = np.fft.fftshift(np.fft.fftfreq(N_samp, d = z_max/N_samp))
t = np.linspace(-pulse_duration*t_max_stdev/2, pulse_duration*t_max_stdev/2, N_samp)
w = np.fft.fftshift(np.fft.fftfreq(N_samp, d = pulse_duration*t_max_stdev/N_samp)) * 2 * const.pi
photon_order = w/w_opt 
g = 49 * 10**-9 * field

### 2.3 Functions modelling the system ###

# def psi_in(t):
#     return np.exp(-t**2/(2*pulse_duration**2)) 

def psi_in(t):
    return np.exp(-t**2/(2*pulse_duration**2)) + np.exp(-(t-pulse_duration)**2/(2*pulse_duration**2)) + np.exp(-(t+pulse_duration)**2/(2*pulse_duration**2)) 

def spec_in_t(t):
    psi_w_squared =np.abs(np.fft.fftshift(np.fft.fft(psi_in(t))))**2
    return psi_w_squared / max(psi_w_squared)

def A(t):
    return np.exp(-2j * np.abs(g) * np.sin(w_opt*t + np.angle(g)))

def psi_t_0(t):
    return A(t) * psi_in(t)

def broadening_fct():
    return np.exp(-w**2/((np.amax(w)-np.amin(w))/1000)**2)

def psi_w_0(t):
    return np.fft.fftshift(np.fft.fft(psi_t_0(t))) 

def spec(t):
    psi_squared = np.abs(psi_w_0(t))**2
    return psi_squared / np.amax(psi_squared)

def propagator(delta_t):
    return np.exp(1j*w*delta_t)

def psi_w_t(t, delta_t):
    # psi_w_0(t) * propagator(delta_t)
    return np.sqrt(spec(t)) * propagator(delta_t)

def spec_t(t, delta_t):
    psi_squared = np.abs(psi_w_t(t, delta_t)**2)
    return psi_squared / np.amax(psi_squared)

def psi_z_t(t, delta_t):
    return np.fft.ifft(psi_w_t(t, delta_t))

def z_vt(delta_t):
    return np.fft.fftshift(np.fft.fftfreq(N_samp, d = np.abs(w[1]-w[0]))) + v_el*delta_t

def electron_density(t, delta_t):
    psi_squared = np.abs(psi_z_t(t, delta_t)**2)
    return psi_squared / np.amax(psi_squared)

### 2.4 Graphing and visualization ###

def drawProgressBar(percent, barLen = 20):
    # percent float from 0 to 1. 
    sys.stdout.write("\r")
    sys.stdout.write("[{:<{}}] {:.0f}%".format("=" * int(barLen * percent), barLen, percent * 100))
    sys.stdout.flush()

def find_peak(spec):
    cutoff = 0.01
    prev_val = 1
    for val in spec:
        if val > cutoff and val > prev_val:
            peak_val = val
            break
        prev_val = val
    index = list(spec).index(peak_val)
    return index

### 2.5 Troubleshooting ###

# fig, ax = plt.subplots()
# ax.plot(t, psi_in(t))
# ax.plot(t, psi_in_corr(t))

fig, ax = plt.subplots()
ax.plt(photon_order, spec(t))

delta_t = 1*10**-10

# fig, ax = plt.subplots()
# ax.plot(z_vt(delta_t), electron_density(t, delta_t))
# # ax.set(ylim=(-1e-20,1e-19))
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(photon_order, spec_t(t, delta_t))
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(w, np.sin(w*10**-16))

plt.show()



