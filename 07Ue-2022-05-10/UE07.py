import numpy as np
import matplotlib.pyplot as plt

# SI units
m_s = 1.989 * 1e30
m_e = 5.972 * 1e24
r_e = 1.4959787 * 1e11
G = 6.67408 * 1e-11
day = 24*3600

# Astro Units
G = G * 1 / r_e**3 * m_e * day **2
m_s = m_s / m_e
m_e = 1
r_e = 1


# Explizites Euler Verfahren

def acc(r):
    norm = np.linalg.norm(r)
    out = - (G * m_s / norm**3 ) * r
    return out
    
v = np.array([0, -0.017326])
r = np.array([1, 0])
dt = 1
pos = [r]


for t in range(int(1 / dt * 365 *2)):
    v_ = v
    v = v + dt * acc(r)
    r = r + dt * v_
    pos.append(r)
pos = np.array(pos)

plt.xlabel('x')
plt.ylabel('y')
plt.title('Earth Orbit integrated with dt = '+str(dt))
plt.plot(pos[:, 0], pos[:, 1])
plt.show()


# Runge-Kutta



v = np.array([0, -0.017326])
r = np.array([1, 0])
dt = 1
pos = [r]


for t in range(int(1 / dt * 365 * 2)):
    k1 = acc(r)
    k1_v = v 

    k2 = acc(r + dt / 2 * k1_v)
    k2_v = v + dt/2*k1

    k3 = acc(r + dt / 2 * k2_v)
    k3_v = v + dt/2*k2

    k4 = acc(r +     dt * k3_v)
    k4_v = v + dt*k3

    phi = 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    phi_v = 1/6 * (k1_v + 2*k2_v + 2*k3_v + k4_v)

    v_ = v
    v = v + dt * phi
    
    r = r + dt * phi_v
    pos.append(r)

pos = np.array(pos)

plt.xlabel('x')
plt.ylabel('y')
plt.title('Earth Orbit integrated with dt = '+str(dt)+' using RK4')
plt.plot(pos[:, 0], pos[:, 1])
plt.show()


#semi-implizierter Euler
v = np.array([0, -0.017326])
r = np.array([1, 0])
dt = 1
pos = [r]

for t in range(int(1 / dt * 365 *2)):
    v = v + dt * acc(r)
    r = r + dt * v
    pos.append(r)
pos = np.array(pos)

plt.xlabel('x')
plt.ylabel('y')
plt.title('Earth Orbit integrated with dt = '+str(dt)+' using semi-implicit Euler')
plt.plot(pos[:, 0], pos[:, 1])
plt.show()