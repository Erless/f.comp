#!/usr/bin/env python
# coding: utf-8

"""
Main for the resolution of the section C in practice 1.
"""

"""
Import of the needed external modules to make the programme work.  
"""
import matplotlib.pyplot as plt
import numpy as np
from Picture_1 import Picture
from Runge_Kutta import RungeKutta_2d # initial form of the object: (to, xo, yo, xf, h, error)


"""
We create empty lists were we are going to store the data that would be plotted.
"""
lista_t = []
lista_x = []
lista_y = []
lista_h = []


"""
The ex_rkf_2d method applies the RKF45 numerical calculation to an object (rk). It calculates the value of the fi and gi.
If the value of rk.aux == True (an assistant optimization controller), it updates the values of position (x) and 
velocity (y) and stores the data in the lists for its representation.
If flag == True (which is set by default, because it increases the precision of the calculations) it also updates 
the value of the integration step (h).
"""
def ex_rkf_2d(rk, flag=(True)):
    rk.fill_lists(lista_t, lista_x, lista_y, lista_h)
    while(rk.t < rk.tf):
        rk.calcular_f0()
        rk.calcular_g0()
        rk.calcular_f1()
        rk.calcular_g1()
        rk.calcular_f2()
        rk.calcular_g2()
        rk.calcular_f3()
        rk.calcular_g3()
        rk.calcular_f4()
        rk.calcular_g4()
        rk.calcular_f5()
        rk.calcular_g5()
        
        rk.aux = True
        
        if (flag == True):
            rk.actualizar_paso_integracion()
            
        if (rk.aux == True):
            rk.actualizar_tiempo()
            rk.actualizar_posicion()
            rk.actualizar_velocidad()
            rk.fill_lists(lista_t, lista_x, lista_y, lista_h)



"""
Defining the functions:
f_1 and f_2 correspond to dy/dx, each one for different problem. They are used in the numerical calculation.
g corresponds to dx/dt and it is the same for both problems.
"""
def f_1(t,x,v): return 4.*(1. - x ** 2.)* v - x
def f_2(t,x,v): return 1. * (1. - x ** 2.) * v - x
def g(t,x,v): return v


"""
We create Two RungeKutta_2d objects with the given conditions, one for solving a different proposed problem.
"""
obj = RungeKutta_2d(0., 2., 0., 50., 0.5, 10 ** (-6))
obj_1 = RungeKutta_2d(0., 0.01, 0.01, 50., 0.5, 10 ** (-6))


"""
Here we give the corresponding function to every object.
"""
obj.f = f_1
obj.g = g

obj_1.f = f_2
obj_1.g = g



"""
We execute the function ex_rkf_2d with the first object.
"""
ex_rkf_2d(obj)


"""
Plotting the results for the first problem:
The first picture contains a representation of the evolution of the x position within time, x(t)
The second picture contains a representation of the evolution of the velocity within time, v(t)
The third picture contains a representation of the phase diagram, v(x)
"""
Picture({'items' : [{'x' : lista_t, 'y' : lista_x, 'legend' : 'x(t)', 'linestyle' : '--'}], 
         'title' : 'X vs t', 'xlabel' : 'time', 'ylabel' : 'x pos.'}).show_plot()
Picture({'items' : [{'x' : lista_t, 'y' : lista_y, 'legend' : 'v(t)', 'linestyle' : '--'}],
         'title' : 'V vs t', 'xlabel' : 'time', 'ylabel' : 'vel.'}).show_plot()
Picture({'items' : [{'x' : lista_x, 'y' : lista_y, 'legend' : 'v(x)', 'linestyle' : '--'}], 
         'title' : 'Phase diagram', 'xlabel' : 'x pos.', 'ylabel' : 'vel.'}).show_plot()


"""
Now we empty the lists, because we use the four lists defined at the begining of the programme to execute the method 
ex_rkf_2d, but they are filled with the first problem data.
"""
lista_t = []
lista_x = []
lista_y = []
lista_h = []



"""
We execute the function ex_rkf_2d with the second object, this time with the update for h set off.
"""
ex_rkf_2d(obj_1)



"""
Plotting the results for the first problem:
The first picture contains a representation of the evolution of the x position within time, x(t)
The second picture contains a representation of the evolution of the velocity within time, v(t)
The third picture contains a representation of the phase diagram, v(x)
"""
Picture({'items' : [{'x' : lista_t, 'y' : lista_x, 'legend' : 'x(t)', 'linestyle' : '--'}], 
         'title' : 'X vs t', 'xlabel' : 'time', 'ylabel' : 'x pos.'}).show_plot()
Picture({'items' : [{'x' : lista_t, 'y' : lista_y, 'legend' : 'y(t)', 'linestyle' : '--'}], 
         'title' : 'V vs t', 'xlabel' : 'time', 'ylabel' : 'vel.'}).show_plot()
Picture({'items' : [{'x' : lista_x, 'y' : lista_y, 'legend' : 'v(x)', 'linestyle' : '--'}], 
         'title' : 'Diagrama de fases', 'xlabel' : 'x pos.', 'ylabel' : 'vel.'}).show_plot()