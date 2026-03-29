import numpy as np
import matplotlib.pyplot as plt

# constants
G = 4 * np.pi**2
M = 1.0

# initial conditions
x0 = 1.0
y0 = 0.0
vx0 = 0.0
vy0 = 2 * np.pi   # circular orbit speed at 1 AU

# time stuff
dt = 0.001
tmax = 3.0
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

# plots




# orbit plot
plt.figure(figsize=(6,6))

# trajectories
plt.plot(x_e, y_e, label='Euler')
plt.plot(x_l, y_l, label='Leapfrog')

# Sun
plt.scatter(0, 0, color='gold', s=200, label='Sun')

# START POINTS (slightly offset so they are visible)
plt.scatter(x_e[0], y_e[0] + 0.02, color='green', s=150, label='Start (offset)')
plt.scatter(x_e[0], y_e[0], color='black', s=40)

# END POINTS
plt.scatter(x_e[-1], y_e[-1], color='red', s=150, label='Euler End')
plt.scatter(x_l[-1], y_l[-1], color='purple', s=150, label='Leapfrog End')

plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.title('Earth Orbit Around the Sun')

plt.axis('equal')
plt.grid()
plt.legend(loc='upper right')

plt.show()




# speed plot
plt.figure()
plt.plot(t, speed_e, label='Euler')
plt.plot(t, speed_l, label='Leapfrog')
plt.xlabel('time [yr]')
plt.ylabel('speed [AU/yr]')
plt.title('Speed vs Time')
plt.grid()
plt.legend()
plt.show()

# energy plot
plt.figure()
plt.plot(t, energy_e, label='Euler')
plt.plot(t, energy_l, label='Leapfrog')
plt.xlabel('time [yr]')
plt.ylabel('specific energy')
plt.title('Energy vs Time')
plt.grid()
plt.legend()
plt.show()




import matplotlib.animation as animation

fig, ax = plt.subplots(figsize=(6,6))

ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')
ax.set_title("Earth Orbit Animation")
ax.grid()

sun, = ax.plot([0], [0], 'yo', markersize=10, label='Sun')

earth_e, = ax.plot([], [], 'bo', label='Euler')
earth_l, = ax.plot([], [], 'ro', label='Leapfrog')

path_e, = ax.plot([], [], 'b-', alpha=0.5)
path_l, = ax.plot([], [], 'r-', alpha=0.5)

ax.legend(loc='upper right')

def update(i):
    earth_e.set_data([x_e[i]], [y_e[i]])
    earth_l.set_data([x_l[i]], [y_l[i]])

    path_e.set_data(x_e[:i+1], y_e[:i+1])
    path_l.set_data(x_l[:i+1], y_l[:i+1])

    return earth_e, earth_l, path_e, path_l

ani = animation.FuncAnimation(
    fig,
    update,
    frames=range(0, n+1, 10),
    interval=30,
    blit=True
)

ani.save("orbit_movie.mp4", writer="ffmpeg", dpi=200)

plt.show()