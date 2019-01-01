#!/usr/bin/env python
# coding: utf-8

"""
Main for the resolution of the section B in practice 1.
"""

"""
Import of the needed external modules to make the programme work.  
"""
import matplotlib.pyplot as plt
import numpy as np
from Picture_1 import Picture
from Runge_Kutta import RungeKutta  # initial form of the object: (xo, yo, xf, h, error)


"""
We create the empty lists were we are going to store the data that would be plotted.
"""
lista_x = []
lista_y = []
lista_h = []



"""
The ex_rkf method applies the RKF45 numerical calculation to an object (rk). It calculates the value of the fi.
If the value of rk.aux == True (an assistant optimization controller), it updates the values of position (x) and 
velocity (y) and stores the data in the lists for its representation.
If flag == True (which is set by default, because it increases the precision of the calculations) it also updates 
the value of the integration step (h).
"""
def ex_rkf(rk, flag=(True)): #rk is the object on which the method is executed
    rk.fill_lists(lista_x, lista_y, lista_h) #Here we add the the first group of values, before starting the calculation

    while(rk.x < rk.xf):
        rk.fill_lists(lista_x, lista_y, lista_h)
        rk.calcular_f0()
        rk.calcular_f1()
        rk.calcular_f2()
        rk.calcular_f3()
        rk.calcular_f4()
        rk.calcular_f5()

        """
        We set the rk.aux to True in order to a proper execution.
        """
        rk.aux = True 
        
        if (flag == True):
            rk.actualizar_paso_integracion() #This method may set rk.aux to False
        if (rk.aux == True):
            rk.actualizar_posicion()
            rk.actualizar_velocidad()
            rk.fill_lists(lista_x, lista_y, lista_h)
        


"""
Defining the functions:
f corresponds to dy/dx. it is used for the numerical calculation.
f_analytic is the problem solving function.
"""
def f(x,y): return (1./(x**2.)+4.*(x-6.)*np.exp(-2.*(x-6.)**2.))*(-1.)
def f_analytic(x): return 1./x + np.exp(-2.*(x-6.)**2.) - np.exp(-50)


"""
We create a RungeKutta object with the given conditions and an array with the x values for the analytical function.
"""
obj = RungeKutta(1, 1, 10, 0.01, 1 * (10 ** -6))
x_analytic = np.arange(1, 10, 0.001)



"""
We execute the method 'ex_rkf' with the obj object in order to find the numerical solution for the given problem
and we run over the x_analytical values, while calculating the analytical solution. 
"""
obj.f = f  #We give obj the desired function
ex_rkf(obj)
y_analytic = f_analytic(x_analytic)

"""
With np.array we fix the format of the lists in order to avoid problems in the calculation of the error
in the numerical solution. Then we calculate it.
"""
x_n = np.array(lista_x)
y_n = np.array(lista_y)
Err=abs(y_n - (1./x_n + np.exp(-2*(x_n - 6)**2)- np.exp(-50)))


"""
Here we create a graphic representation for each given question for a comfortable analysis of the results.
The first Picture is a comparison between the numeric and the analytic results (y(x)).
The second Picture shows a representation of the evolution of the integration step with x (h(x)) and the last Picture
shows the error for the numerical calculation.
"""
Picture({'items' : [{'x' : lista_x, 'y' : lista_y, 'legend' : 'RKF-4/5 sol. ', 'linestyle' : '--'},
                   {'x' : x_analytic, 'y' : y_analytic, 'legend' : 'Analytical sol.'}],
                     'title' : 'y vs x', 'xlabel':'pos x', 'ylabel':'pos y'}).show_plot()


Picture({'items' : [{'x' : lista_x, 'y' : lista_h, 'legend' : 'h'}],
         'title' : 'h vs x', 'xlabel':'pos x', 'ylabel':'h'}).show_plot()


Picture({'items' : [{'x' : lista_x, 'y' : Err, 'legend' : 'Error'}],
         'title' : 'Error vs x', 'xlabel':'pos x', 'ylabel':'Err'}).show_plot()