#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW10 - Problem 5. Multiple fission algae model, using optimization

Created on Sun Feb 14 00:40:51 2021
@author: eduardo
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def logistic_model(K, P0, r, t):
    """
    Logistic model for population growth

    Parameters
    ----------
    K : float
        carrying capacity.
    P0 : int
        initial population.
    r : float
        low density growth rate.
    t : array
        actual time vector.

    Returns
    -------
    array
        population for each time step.

    """
    return (K * P0 * np.exp(r*t))/(K + P0 * (np.exp(r*t) - 1))
    
def sum_squared_error(x, t, P):
    """
    Sum of squared errors

    Parameters
    ----------
    x : array
        array containing  parameters K, P0, and r (in that order).
    t : array
        actual time vector.
    P : array
        population to compare against.

    Returns
    -------
    SSE : float
        sum of squared errors.

    """
    assert len(t) == len (P), "Length of t and P are not equal"
    K, P0, r = x[0], x[1], x[2]  # unzip the parameters
    Pop = logistic_model(K, P0, r, t)  # calculate population with logistic model (analytical)
    SS = pow(Pop - P, 2)  # squared errors between analytical and numerical populations
    SSE = sum(SS)
    return SSE

kN = 3e-5  # Nitrogen, mole/L
kP = 1e-5  # Phophorous, mole/L
kCO2 = 2.5e-5  # Carbon dioxide, mole/L
time_step = 1 # time step, hours
C0 = 0.05  # Initial concentration, mg/L
max_light = 800  # W/m2
k = 0.002  # Light mu multiplier, W/m2
mu_night = -0.1  # 1/day
number_of_days = 20  # days
max_conc = 2 #  mg/L
N0 = 0.1  # Nitrogen initial, moles/L
P0 = 0.01  # Phophorous initial, moles/L
CO20 = 0.1  # Carbon dioxide initial, moles/L

hpd = 24  # hours per day
N_moles = 16
P_moles = 1
CO2_moles = 124

h = time_step * (1/hpd)  # convert time step into days
time_steps = int(number_of_days / h)  # get time steps in the number of days
# This is grams of algae per mole of P, or 16 moles of N
algae_conc = 12*106 + 1*263 + 16*110 + 14*16 + 31  # total g per mole

time = np.linspace(0, (time_steps*h), time_steps+1)  # time vector

C = np.ones(time_steps+1) * C0
N = np.ones(time_steps+1) * N0
P = np.ones(time_steps+1) * P0
CO2 = np.ones(time_steps+1) * CO20
I0 = np.zeros(time_steps+1)
mu = np.ones(time_steps+1) * mu_night

for i in range(1, time_steps+1):
    I0[i] = max_light * np.sin(2 * np.pi * (i-6)/hpd)
    I0_effective = I0[i] * (max_conc - C[i-1]) / max_conc
    
    mu[i] = I0_effective * k * N[i-1]/(kN + N[i-1]) * P[i-1]/(kP + P[i-1]) * CO2[i-1]/(kCO2 + CO2[i-1])
    
    mu[i] = mu_night if mu[i] < 0 else mu[i]
    
    deltaC = C[i-1] * mu[i] * h
    C[i] = C[i-1] + deltaC
    
    if mu[i] > 0:
        N[i] = N[i-1] - deltaC / algae_conc * N_moles
        P[i] = P[i-1] - deltaC / algae_conc * P_moles
        CO2[i] = CO2[i-1] - deltaC / algae_conc * CO2_moles
    else:
        N[i] = N[i-1]
        P[i] = P[i-1]
        CO2[i] = CO2[i-1]

# Use 'minimize' to find the coefficients K, P0, r
x0 = [0.1, 0.1, 0.1]  # initial guess for K, P0, r
res = minimize(sum_squared_error, x0, args=(time, C))
K, P0, r = res.x[0], res.x[1], res.x[2]  # get the results
print('K=', K, 'P0=', P0, 'r=', r)

# Make predictions using the analytical model and numpy array
# coefficients were obtained from Excel using Solver
Cexp = logistic_model(K, C0, r, time)

growth_script = pd.DataFrame(columns = ['Time', 'I0', 'mu', 'Conc', 'N', 'P', 'CO2'])
growth_script['Time'] = time
growth_script['I0'] = I0
growth_script['mu'] = mu
growth_script['Conc'] = C
growth_script['N'] = N
growth_script['P'] = P
growth_script['CO2'] = CO2

plt.figure(0)
plt.plot(growth_script['Time'], growth_script['Conc'], 'k-', label='Conc')
plt.plot(growth_script['Time'], Cexp, 'r-', label='Conc (exp)')
# plt.plot(growth_script['Time'], growth_script['N'], 'b--', label='N')
# plt.plot(growth_script['Time'], growth_script['P'], 'r-.', label='P')
# plt.plot(growth_script['Time'], growth_script['CO2'], 'm:', label='$CO_2$')
# plt.plot(growth_script['Time'], growth_script['mu'], 'y-', label='$\mu$')
plt.legend(loc='best', fancybox=True)
plt.xlabel('Time [days]')
plt.ylabel('Concentration $K,P,CO_2$ [mole/L], C [mg/L]')
# plt.savefig('p5_algae_model_%ddays.png' % number_of_days, dpi=300, bbox_inches='tight')
plt.show()
