import numpy as np
import matplotlib.pyplot as plt
import time

class Param:
    def __init__(self, *args):
        pass

def nOrd(WL):
    WL  = WL/1e3
    return np.sqrt(2.7359 + 0.01878 /(WL**2 - 0.01822) - 0.01354*WL**2)

def nExtra(WL):
    WL = WL/1e3
    return np.sqrt(2.3753 + 0.01224 /(WL**2 - 0.01667) - 0.01516*WL**2)

def nExtraTeta(WL, TETA):
    nO = nOrd(WL)
    nE = nExtra(WL)
    return nO * nE/np.sqrt(nE**2 * cosd(TETA)**2 + nO**2 * sind(TETA)**2)

def sind(angle):
    return(np.sin(angle / 180 * np.pi))

def cosd(angle):
    return(np.cos(angle / 180 * np.pi))

Tstart = time.time()

Par = Param()
Par.N1, Par.N2 = int(1e3), int(1e3)
tetaCr = 48 + 0.5
L = 2e-3

WL = np.linspace(600, 800, Par.N1) # Диапазон длин волн
TETA = np.linspace(-5, 5, Par.N2) # Диапазон углов

WLs, TETAs = np.meshgrid(WL, TETA) #

c = 3e8
pump = 355e0

wp = c / pump * 1e9
ws = c / WLs * 1e9
wi = wp-ws

n_idle = nOrd(c/wi*1e9)
n_pump = nExtraTeta(pump, tetaCr)
n_signal_1 = nExtraTeta(c/ws*1e9, tetaCr-TETAs)
n_signal_2 = nExtraTeta(c/ws*1e9, tetaCr+TETAs)

delta_kp = lambda: n_pump*wp/c
delta_ks = lambda n1: n1*ws/c * np.sqrt( 1- ( sind(TETAs)/n1 )**2 )
delta_ki = lambda n2: n2*(wp-ws)/c * np.sqrt( 1 - ( ws*sind(TETAs)/n2/(wp-ws) )**2 )
delta_K = lambda n1, n2: delta_kp() - delta_ki(n1) - delta_ks(n2)

F_o = L * np.sinc( L/2 * delta_K(n_idle, n_signal_1) )**2
F_e = L * np.sinc( L/2 * delta_K(n_signal_2, n_idle) )**2

SPDC_spectre = 1e3 * np.abs(F_e + F_o)

Tend = time.time()

fig ,ax = plt.subplots()
pc = ax.pcolormesh(WL, TETA, SPDC_spectre)
plt.colorbar(pc)
plt.show()

print(f'{Tend - Tstart} seconds')
