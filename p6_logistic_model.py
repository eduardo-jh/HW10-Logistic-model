# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW10 - Problem 6. Logistic model

Created on Sat Feb 13 22:56:55 2021
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt

dt = 1/24  # time step, in days
K = 3000  # maximum carrying capacity
r = 0.4  # low density growth rate
P0 = 10  # initial population
steps = 25

t = np.linspace(0, steps, int(steps/dt)+1)  # time vector
print(t)

# initialize the population vectors
P = np.zeros(len(t))
Pana = np.zeros(len(t))
P[0], Pana[0] = P0, P0

for i in range(1, len(t)):
    P[i] = P[i-1] + r * P[i-1] * (1 - P[i-1]/K)*dt
    Pana[i] = (K * P0 * np.exp(r*t[i]))/(K + P0 * (np.exp(r*t[i]) - 1))

# Make predictions using the analytical model and numpy array
Pexp = (K * P0 * np.exp(r*t))/(K + P0 * (np.exp(r*t) - 1))

plt.plot(t, P, 'bx', t, Pana, 'g+', t, Pexp, 'r-')
plt.legend(['Euler (dt=%.1f h)' % (dt*24), 'Analytic (loop)', 'Analytic (eq.)'], loc='best')
plt.xlabel('Time (days)')
plt.ylabel('Population')
plt.savefig('p6_logistic_model_%ddays.png' % steps, dpi=300, bbox_inches='tight')
plt.show()
