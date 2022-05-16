
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Units

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


# Physics

def acc(r):
    norm = np.linalg.norm(r)
    out = - (G * m_s / norm**3 ) * r
    return out


# Explizites Euler Verfahren

def exp_Eu(r0, v0, dt, T):
    r = r0
    v = v0
    pos = [r0]
    vel = [v0]
    for t in range(int(1 / dt * T)):
        v_ = v
        v = v + dt * acc(r)
        r = r + dt * v_
        pos.append(r)
        vel.append(v)

    pos = np.array(pos)
    v = np.array(v)
    return pos, v

r0 = np.array([1, 0])
v0 = np.array([0, -2*np.pi / 365])
dt1 = 1
dt2 = 0.01
T = 365

pos1 = exp_Eu(r0, v0, dt1, T)[0]
pos2 = exp_Eu(r0, v0, dt2, T)[0]

plt.xlabel('x')
plt.ylabel('y')
plt.title('Earth Orbit integrated using the explicit euler method')
plt.plot(pos1[:, 0], pos1[:, 1], label="dt = "+str(dt1)+' days')
plt.plot(pos2[:, 0], pos2[:, 1], label="dt = "+str(dt2)+' days')
plt.legend()
plt.show()


# Runge-Kutta

def RK4(r0, v0, dt, T):
    r = r0
    v = v0
    pos = [r0]
    vel = [v0]

    for t in range(int(1 / dt * T)):
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
        
        v = v + dt * phi
        r = r + dt * phi_v

        pos.append(r)
        vel.append(v)

    pos = np.array(pos)
    vel = np.array(vel)
    return pos, vel

pos1 = RK4(r0, v0, dt1, T)[0]
pos2 = RK4(r0, v0, dt2, T)[0]

plt.xlabel('x')
plt.ylabel('y')
plt.title('Earth Orbit integrated using the Runge-Kutta method')
plt.plot(pos1[:, 0], pos1[:, 1], label="dt = "+str(dt1)+' days')
plt.plot(pos2[:, 0], pos2[:, 1], label="dt = "+str(dt2)+' days')
plt.legend()
plt.show()


#semi-implizierter Euler

def SIM_Eu(r0, v0, dt, T):
    r = r0
    v = v0
    pos = [r0]
    vel = [v0]
    for t in range(int(1 / dt * T)):
        v = v + dt * acc(r)
        r = r + dt * v
        pos.append(r)
        vel.append(v)
    pos = np.array(pos)
    vel = np.array(vel)
    return pos, vel

v0 = np.array([0, -0.017326])
r0 = np.array([1, 0])
dt = 1
T = 365

pos = SIM_Eu(r0, v0, dt, T)[0]

plt.xlabel('x')
plt.ylabel('y')
plt.title('Earth Orbit integrated with dt = '+str(dt)+' using semi-implicit Euler')
plt.plot(pos[:, 0], pos[:, 1])
plt.show()


#impliziertes Eulerverfahren 

def imp_Eu(r0, v0, dt, T):
    r = r0
    v = v0
    pos = [r0]
    vel = [v0]
    for t in range(int(1 / dt * T)):
        r = r + dt * v
        v = v + dt * acc(r)
        pos.append(r)
        vel.append(v)
    pos = np.array(pos)
    vel = np.array(vel)
    return pos, vel

r0 = np.array([1, 0])
v0 = np.array([0, -2*np.pi / 365])
dt1 = 1
dt2 = 0.01
T = 365

pos1 = imp_Eu(r0, v0, dt1, T)[0]
pos2 = imp_Eu(r0, v0, dt2, T)[0]

plt.xlabel('x')
plt.ylabel('y')
plt.title('Earth Orbit integrated using the implicit euler method')
plt.plot(pos1[:, 0], pos1[:, 1], label="dt = "+str(dt1)+' days')
plt.plot(pos2[:, 0], pos2[:, 1], label="dt = "+str(dt2)+' days')
plt.legend()
plt.show()


# Energien 

def E(pos, vel):
    E = m_e*vel/2 - G*m_e*m_s/np.abs(pos)
    return E

dt = 0.01
r0 = np.array([1, 0])
v0 = np.array([0, -2*np.pi / 365])


pos_impE = imp_Eu(r0, v0, dt1, T)[0]
vel_impE = imp_Eu(r0, v0, dt1, T)[1]

pos_expE = exp_Eu(r0, v0, dt1, T)[0]
vel_expE = exp_Eu(r0, v0, dt2, T)[1]

pos_RK4 = RK4(r0, v0, dt1, T)[0]
vel_RK4 = RK4(r0, v0, dt2, T)[1]


E_impE=E(pos_impE, vel_impE)
E_expE=E(pos_expE, vel_expE)
#E_RK4=E(pos_RK4, vel_RK4)
print(vel_RK4.shape)
print(pos_RK4.shape)
print(vel_RK4)


plt.xlabel('t-steps')
plt.ylabel('E')
plt.title('Energies')
plt.plot(E_impE, label='implicit Euler')
plt.plot(E_expE, label='explicit Euler')
plt.legend()
plt.show()


# Transfer
# We want to find a velocity, such that the maximum radial distance is R_m = 1.5
# according to kepler the orbit should have a period about 1.837 < 2 times the period of
# the earth orbit

R_mars = 1.5

r0 = np.array([1, 0])
v0 = np.array([0, -2*np.pi / 365])
dt1 = 0.1
dv = v0[1] / 100

pos = RK4(r0, v0, dt1, 2*T)[0]
r_max = np.max(np.abs(pos))
radii = []

while r_max < R_mars:
    v0[1] = v0[1] + dv
    pos = RK4(r0, v0, dt1, 2*T)[0]
    r_max = np.max(np.abs(pos))
    #print(r_max, v0)
    radii.append([v0[1], r_max])

radii = np.array(radii)

fit = interp1d(radii[:, 1], radii[:, 0], kind='cubic')
V_mars = fit(R_mars)

pos_mars = RK4(r0, np.array([0, V_mars]), dt1, 2*T)[0]
t_max = np.where(np.abs(pos_mars) == np.max(np.abs(pos_mars)))

print('Tage bis zum Mars: ', t_max[0] * dt1)
