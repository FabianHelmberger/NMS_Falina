#/usr/bin/python3
#-*- coding: utf-8 -*-

'''
Ue3
---------------

This module is meant to provide solutions to the questions asked in exercise
seven.
'''

import numpy as np
from functools import wraps
import matplotlib.pyplot as plt

def accelerate(r):
    acc = -GAU*r*SUNEARTH/np.linalg.norm(r)**3
    return acc

def step_iterator(n=365):
    def wrap(function):
        def wrapped_method(r,v,h):
            position = []
            velocity = []
            iteration = 0
            while True:
                position.append(r)
                velocity.append(v)
                r, v = function(r, v, h)
                iteration += 1
                if iteration == n/h:
                    return position, velocity
        return wrapped_method
    return wrap

@step_iterator()
def explicit_euler(r, v, h):
    _r = r+v*h
    _v = v+accelerate(r)*h
    return _r, _v
           
@step_iterator()
def rk4(r,v,h):
    k1 = v 
    l1 = accelerate(r)
    k2 = v + h/2*l1
    l2 = accelerate(r+h/2*k1)
    k3 = v + h/2*l2
    l3 = accelerate(r+h/2*k2)
    k4 = v+h*l3
    l4 = accelerate(r+h*k3)
    phi_v = 1/6*(k1 + 2*k2 + 2*k3 + k4)
    phi = 1/6*(l1 + 2*l2 + 2*l3 + l4)
    v = v + phi*h
    r = r + phi_v*h
    return r,v

@step_iterator()
def semi_euler(r, v, h):
    r = r+v*h
    v = v+accelerate(r)*h
    return r, v

@step_iterator(700)
def mars_semi_euler(r, v, h):
    r = r+v*h
    v = v+accelerate(r)*h
    return r, v

def plot_solutions(**kwargs):
    for key, value in kwargs.items():
        if 'title' in key:
            plt.title(value)
        elif 'dt' in key:
            x = list(zip(*value))[0]
            y = list(zip(*value))[1]
            plt.plot(x,y, label = key)
        else:
            plt.plot(value, label = key)
    plt.legend()
    plt.show()
    
def energy(positions, velocities):
    energy = []
    for r, v in zip(positions, velocities):
        e = np.linalg.norm(v)**2/2 - GAU*SUNEARTH/np.linalg.norm(r)
        energy.append(e)
    return energy

if __name__ == '__main__':

    MS = 1.989e30
    ME = 5.972e24
    RE = 1.4959787e11
    G = 6.67408e-11
    SUNEARTH = MS/ME
    GAU = G*1/RE**3*ME*86400**2
    
    r = np.array([1,0])
    v = np.array([0, -0.017326])
    position1, velocity1 = explicit_euler(r,v,1)
    position2, velocity2 = explicit_euler(r,v,0.01)
    explicit_euler_energies = energy(position2, velocity2)
    plot_solutions(dt_86400_seconds= position1, dt_8640_seconds=position2, 
        title = 'Annual orbit of the Earth approximated using Explicit Euler method')

    position1, velocity2 = rk4(r,v,1)
    position2, velocity2 = rk4(r,v,0.01)
    rk4_energies = energy(position2, velocity2)
    plot_solutions(dt_86400_seconds= position1, dt_8640_seconds=position2, 
        title = 'Annual orbit of the Earth approximated using RK4 method')

    position1, velocity2 = semi_euler(r,v,1)
    position2, velocity2 = semi_euler(r,v,0.01)
    implicit_euler_energies = energy(position2, velocity2)
    plot_solutions(dt_86400_seconds= position1, dt_8640_seconds=position2, 
        title = 'Annual orbit of the Earth approximated using Implicit Euler method')

    plot_solutions(energies_explicit_euler_method = explicit_euler_energies,
        runge_kutta_energies = rk4_energies, 
        energies_implicit_euler = implicit_euler_energies,
        title = 'Energy of the entire System vs time calculated with'\
            ' different approximations using dt=0.01 day')

    positions, velocity = mars_semi_euler(np.array([1,0]), np.array([0,-0.02]), 0.01)
    plot_solutions(dt_8640_seconds_Mars_orbit = positions, 
        dt_8640_seconds_Earth_orbit = position2,
        title = 'No title added yet')

    
