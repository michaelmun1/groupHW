# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 21:38:25 2026

@author: micha
"""

import numpy as np
import matplotlib.pyplot as plt

# constants
G = 4 * np.pi**2
M = 1.0

# initial conditions
x0 = 1.0
y0 = 0.0
vx0 = 0.0
vy0 = 0.8 * 2 * np.pi   # 0.8 times circular orbit speed

# time stuff
dt = 0.001
tmax = 5.0
n = int(tmax / dt)

t = np.zeros(n+1)

#  EULER 
x_e = np.zeros(n+1)
y_e = np.zeros(n+1)
vx_e = np.zeros(n+1)
vy_e = np.zeros(n+1)
speed_e = np.zeros(n+1)
energy_e = np.zeros(n+1)

x_e[0] = x0
y_e[0] = y0
vx_e[0] = vx0
vy_e[0] = vy0

for i in range(n):
    r = np.sqrt(x_e[i]**2 + y_e[i]**2)

    ax = -G * M * x_e[i] / r**3
    ay = -G * M * y_e[i] / r**3

    vx_e[i+1] = vx_e[i] + ax * dt
    vy_e[i+1] = vy_e[i] + ay * dt

    x_e[i+1] = x_e[i] + vx_e[i] * dt
    y_e[i+1] = y_e[i] + vy_e[i] * dt

    t[i+1] = t[i] + dt

for i in range(n+1):
    r = np.sqrt(x_e[i]**2 + y_e[i]**2)
    v = np.sqrt(vx_e[i]**2 + vy_e[i]**2)
    speed_e[i] = v
    energy_e[i] = 0.5 * v**2 - G * M / r

#  LEAPFROG 
x_l = np.zeros(n+1)
y_l = np.zeros(n+1)
vx_l = np.zeros(n+1)
vy_l = np.zeros(n+1)
speed_l = np.zeros(n+1)
energy_l = np.zeros(n+1)

x_l[0] = x0
y_l[0] = y0
vx_l[0] = vx0
vy_l[0] = vy0

for i in range(n):
    r = np.sqrt(x_l[i]**2 + y_l[i]**2)

    ax = -G * M * x_l[i] / r**3
    ay = -G * M * y_l[i] / r**3

    vx_half = vx_l[i] + 0.5 * ax * dt
    vy_half = vy_l[i] + 0.5 * ay * dt

    x_l[i+1] = x_l[i] + vx_half * dt
    y_l[i+1] = y_l[i] + vy_half * dt

    r_new = np.sqrt(x_l[i+1]**2 + y_l[i+1]**2)
    ax_new = -G * M * x_l[i+1] / r_new**3
    ay_new = -G * M * y_l[i+1] / r_new**3

    vx_l[i+1] = vx_half + 0.5 * ax_new * dt
    vy_l[i+1] = vy_half + 0.5 * ay_new * dt

for i in range(n+1):
    r = np.sqrt(x_l[i]**2 + y_l[i]**2)
    v = np.sqrt(vx_l[i]**2 + vy_l[i]**2)
    speed_l[i] = v
    energy_l[i] = 0.5 * v**2 - G * M / r





#  TRAJECTORY PLOT 
plt.figure(figsize=(6,6))

plt.plot(x_e, y_e, label='Euler')
plt.plot(x_l, y_l, label='Leapfrog')

plt.scatter(0, 0, color='gold', s=200, label='Sun')
plt.scatter(x_e[0], y_e[0] + 0.03, color='green', s=150, label='Start')
plt.scatter(x_e[0], y_e[0], color='black', s=40)

plt.scatter(x_e[-1], y_e[-1], color='red', s=150, label='Euler End')
plt.scatter(x_l[-1], y_l[-1], color='purple', s=150, label='Leapfrog End')

plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.title('Q2 Trajectory: v0 = 0.8 v_orb')
plt.axis('equal')
plt.grid()
plt.legend(loc='upper right')
plt.savefig('q2_trajectory.png', dpi=300)
plt.show()






#  EULER ENERGY 
plt.figure()
plt.plot(t, energy_e)
plt.xlabel('time [yr]')
plt.ylabel('specific energy')
plt.title('Q2 Euler Energy vs Time')
plt.grid()
plt.savefig('q2_energy_euler.png', dpi=300)
plt.show()

#  LEAPFROG ENERGY 
plt.figure()
plt.plot(t, energy_l)
plt.xlabel('time [yr]')
plt.ylabel('specific energy')
plt.title('Q2 Leapfrog Energy vs Time')
plt.grid()
plt.savefig('q2_energy_leapfrog.png', dpi=300)
plt.show()













