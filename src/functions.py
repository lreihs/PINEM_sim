import numpy as np
from read_input import *

def psi0():
    return (2*np.pi*sigma_z**2)**-0.25

def psi_z_0():
    return psi0()*np.exp(-z**2/(4*sigma_z**2))

def psiD_z(tD):
    LD = v0*tD
    disp_param = hbar/(2*m_e*lorentz_factor**3*sigma_z**2*v0)
    aD = 1+1j*disp_param*LD
    return psi0()*aD**-0.5*np.exp(-(z-LD)**2/(4*sigma_z**2*aD))

def psi1_z(tD):
    LD = v0*tD
    return psiD_z(tD)*np.exp(1j*np.abs(g)*np.sin(k_p*(z-LD)+np.angle(g)))

def psi1_k(tD):
    return np.fft.fftshift(np.fft.fft(psi1_z(tD)))

def spec1(tD):
    return np.abs(psi1_k(tD))**2

def psiL_k(tD, L):
    E=hbar*k*v0
    return psi1_k(tD)*np.exp(-1j*np.pi*E**2/(hbar*w_opt)**2*L/LQR)

def psiL_z(tD, L):
    return np.fft.ifft(psiL_k(tD, L))

def electron_density(tD, L):
    return np.abs(psiL_z(tD, L))**2




