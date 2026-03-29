# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 21:45:56 2026

@author: micha
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# constants
G = 4 * np.pi**2
M_sun = 1.0
M_jup = 0.0009543

# conversion
AUyr_to_kms = 4.743

#  INITIAL SETUP 

# Jupiter orbit
rJ = 5.2
vJ = np.sqrt(G * M_sun / rJ)



# choose Jupiter's starting angle - gotta djust this where needed
phiJ0_deg = 80
phiJ0 = np.radians(phiJ0_deg)

xj0 = rJ * np.cos(phiJ0)
yj0 = rJ * np.sin(phiJ0)

vxj0 = -vJ * np.sin(phiJ0)
vyj0 =  vJ * np.cos(phiJ0)

# Voyager starts at Earth's orbit
xv0 = 1.0
yv0 = 0.0

# escape speed at Earth orbit
vesc = np.sqrt(2 * G * M_sun / 1.0)

# launch angle measured from straight up (+y direction)
# 0 means purely tangential
launch_angle_deg = 8
launch_angle = np.radians(launch_angle_deg)

# Voyager initial velocity
vxv0 = vesc * np.sin(launch_angle)
vyv0 = vesc * np.cos(launch_angle)

#  TIME 

dt = 0.001
tmax = 3.0
n = int(tmax / dt)

t = np.zeros(n + 1)

# ARRAYS 

# Jupiter
xj = np.zeros(n + 1)
yj = np.zeros(n + 1)
vxj = np.zeros(n + 1)
vyj = np.zeros(n + 1)

# Voyager
xv = np.zeros(n + 1)
yv = np.zeros(n + 1)
vxv = np.zeros(n + 1)
vyv = np.zeros(n + 1)

# speed and distance tracking
voy_speed = np.zeros(n + 1)
sep = np.zeros(n + 1)

# initial values
xj[0] = xj0
yj[0] = yj0
vxj[0] = vxj0
vyj[0] = vyj0

xv[0] = xv0
yv[0] = yv0
vxv[0] = vxv0
vyv[0] = vyv0

voy_speed[0] = np.sqrt(vxv0**2 + vyv0**2)
sep[0] = np.sqrt((xv0 - xj0)**2 + (yv0 - yj0)**2)

# ---------------- LEAPFROG LOOP ----------------

for i in range(n):
    # Jupiter acceleration from Sun only
    rj = np.sqrt(xj[i]**2 + yj[i]**2)
    axj = -G * M_sun * xj[i] / rj**3
    ayj = -G * M_sun * yj[i] / rj**3

    # Voyager acceleration from Sun
    rv = np.sqrt(xv[i]**2 + yv[i]**2)
    axv_sun = -G * M_sun * xv[i] / rv**3
    ayv_sun = -G * M_sun * yv[i] / rv**3

    # Voyager acceleration from Jupiter
    dx = xv[i] - xj[i]
    dy = yv[i] - yj[i]
    r_vj = np.sqrt(dx**2 + dy**2)

    axv_jup = -G * M_jup * dx / r_vj**3
    ayv_jup = -G * M_jup * dy / r_vj**3

    axv = axv_sun + axv_jup
    ayv = ayv_sun + ayv_jup

    # half-step velocities
    vxj_half = vxj[i] + 0.5 * axj * dt
    vyj_half = vyj[i] + 0.5 * ayj * dt

    vxv_half = vxv[i] + 0.5 * axv * dt
    vyv_half = vyv[i] + 0.5 * ayv * dt

    # position update
    xj[i+1] = xj[i] + vxj_half * dt
    yj[i+1] = yj[i] + vyj_half * dt

    xv[i+1] = xv[i] + vxv_half * dt
    yv[i+1] = yv[i] + vyv_half * dt

    # new Jupiter acceleration
    rj_new = np.sqrt(xj[i+1]**2 + yj[i+1]**2)
    axj_new = -G * M_sun * xj[i+1] / rj_new**3
    ayj_new = -G * M_sun * yj[i+1] / rj_new**3

    # new Voyager acceleration
    rv_new = np.sqrt(xv[i+1]**2 + yv[i+1]**2)
    axv_sun_new = -G * M_sun * xv[i+1] / rv_new**3
    ayv_sun_new = -G * M_sun * yv[i+1] / rv_new**3

    dx_new = xv[i+1] - xj[i+1]
    dy_new = yv[i+1] - yj[i+1]
    r_vj_new = np.sqrt(dx_new**2 + dy_new**2)

    axv_jup_new = -G * M_jup * dx_new / r_vj_new**3
    ayv_jup_new = -G * M_jup * dy_new / r_vj_new**3

    axv_new = axv_sun_new + axv_jup_new
    ayv_new = ayv_sun_new + ayv_jup_new

    # full-step velocities
    vxj[i+1] = vxj_half + 0.5 * axj_new * dt
    vyj[i+1] = vyj_half + 0.5 * ayj_new * dt

    vxv[i+1] = vxv_half + 0.5 * axv_new * dt
    vyv[i+1] = vyv_half + 0.5 * ayv_new * dt

    t[i+1] = t[i] + dt

    voy_speed[i+1] = np.sqrt(vxv[i+1]**2 + vyv[i+1]**2)
    sep[i+1] = np.sqrt((xv[i+1] - xj[i+1])**2 + (yv[i+1] - yj[i+1])**2)
    
    
    

# CLOSEST APPROACH

imin = np.argmin(sep)
closest_distance_AU = sep[imin]
closest_time_yr = t[imin]
closest_time_months = 12 * closest_time_yr

print("Closest approach distance =", closest_distance_AU, "AU")
print("Closest approach time =", closest_time_yr, "yr")
print("Closest approach time =", closest_time_months, "months")



#  TRAJECTORY PLOT 

plt.figure(figsize=(7,7))
plt.plot(xj, yj, label='Jupiter')
plt.plot(xv, yv, label='Voyager 2')

plt.scatter(0, 0, color='gold', s=200, label='Sun')
plt.scatter(xj[0], yj[0], color='orange', s=100, label='Jupiter Start')
plt.scatter(xv[0], yv[0], color='green', s=100, label='Voyager Start')

plt.scatter(xj[imin], yj[imin], color='red', s=120, label='Jupiter at Closest Approach')
plt.scatter(xv[imin], yv[imin], color='purple', s=120, label='Voyager at Closest Approach')

plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.title('Voyager 2 Gravity Assist with Jupiter')
plt.axis('equal')
plt.grid()
plt.legend(loc='upper right')
plt.savefig('q3_trajectory.png', dpi=300)
plt.show()

# ---------------- SPEED PLOT ----------------

plt.figure()
plt.plot(t * 12, voy_speed * AUyr_to_kms)
plt.xlabel('time [months]')
plt.ylabel('Voyager speed [km/s]')
plt.title('Voyager Speed vs Time')
plt.grid()
plt.savefig('q3_speed.png', dpi=300)
plt.show()

# ---------------- ANIMATION ----------------

fig, ax = plt.subplots(figsize=(7,7))
ax.set_xlim(-8, 8)
ax.set_ylim(-8, 8)
ax.set_aspect('equal')
ax.grid()
ax.set_title('Voyager 2 and Jupiter')

sun, = ax.plot([0], [0], 'yo', markersize=10)

jup_dot, = ax.plot([], [], 'o', color='orange', markersize=8, label='Jupiter')
voy_dot, = ax.plot([], [], 'o', color='blue', markersize=6, label='Voyager 2')

jup_path, = ax.plot([], [], color='orange', alpha=0.5)
voy_path, = ax.plot([], [], color='blue', alpha=0.5)

text_time = ax.text(-7.5, 7.2, '', fontsize=10)
text_speed = ax.text(-7.5, 6.5, '', fontsize=10)

ax.legend(loc='upper right')

def update(i):
    jup_dot.set_data([xj[i]], [yj[i]])
    voy_dot.set_data([xv[i]], [yv[i]])

    jup_path.set_data(xj[:i+1], yj[:i+1])
    voy_path.set_data(xv[:i+1], yv[:i+1])

    text_time.set_text('Time = %.1f months' % (12 * t[i]))
    text_speed.set_text('Voyager speed = %.2f km/s' % (voy_speed[i] * AUyr_to_kms))

    return jup_dot, voy_dot, jup_path, voy_path, text_time, text_speed

ani = animation.FuncAnimation(
    fig,
    update,
    frames=range(0, n + 1, 10),
    interval=30,
    blit=True
)

# try mp4 first
ani.save('q3_voyager_movie.mp4', writer='ffmpeg', dpi=200)

plt.show()




closest_distance_km = closest_distance_AU * 1.496e8

print("Closest approach distance =", closest_distance_AU, "AU")
print("Closest approach distance =", closest_distance_km, "km")
print("Closest approach time =", closest_time_yr, "yr")
print("Closest approach time =", closest_time_months, "months")