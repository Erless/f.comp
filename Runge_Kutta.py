#!/usr/bin/env python
# coding: utf-8

"""
Methods RungeKutta and RungeKutta_2d.
They allow to create objects that apply the Runge Kutta 4/5 numerical calculation, with or without the adaptive
step added by Fehlberg in order to find the solution to 1st order ODEs in 1-D  or 2-D respectively.
"""



class RungeKutta(object):
    """"
    Builder of RungeKutta class. We give the needed data in order to solve the problem:
    x: Initial x position. It is updated when the calculations are executed.
    y: Initial velocity.  It is updated when the calculations are executed.
    xf: Final x position. Given in order to stop the calculation.
    h: Integration step. It may be updated when the calculations are executed if
    step update method (actualizar_paso_integracion) is active.
    error: Desired error.
    """
    def __init__(self, x, y, xf, h, error, *args, **kwargs):
        self.x = x
        self.y = y
        self.xf = xf
        self.h = h
        self.error = error
        self.f = None
        self.aux = True  #Auxiliar used to optimize the calculations.


    """
    In the next methods, calcular_fi, we implement the calculation of the intermediate functions needed for the 
    calculation of yn+1 and the update of the integration step. 
    """
    
    def calcular_f0(self):
        self.f0 = self.f(self.x, self.y)
    
    def calcular_f1(self):    
        self.f1 = self.f(self.x + self.h / 4., self.y + self.h / 4. * self.f0)
    
    def calcular_f2(self):
        self.f2 = self.f(self.x + 3. * self.h / 8., self.y + 3. * self.h / 32 *
                         self.f0 + 9. * self.h / 32. * self.f1)
        
    def calcular_f3(self):
        self.f3 = self.f(self.x + 12. * self.h/13.2, self.y + 1932. *
                         self.h / 2197. * self.f0 - 7200. * self.h/2197. *
                         self.f1 + 7296. * self.h / 2197. * self.f2)
        
    def calcular_f4(self):
        self.f4 =self.f(self.x + self.h, self.y + 439. * self.h / 216. * 
                        self.f0 - 8. * self.h * self.f1 + 3680. * self.h / 513. *
                        self.f2 - 845. * self.h / 4104. * self.f3)
        
    def calcular_f5(self):
        self.f5 = self.f(self.x + self.h / 2., self.y - 8. * self.h / 27. * self.f0 + 2. *
                         self.h * self.f1 - 3544. * self.h / 2565. * self.f2 + 1859. * self.h / 4104. *
                         self.f3 - 11. * self.h / 40. * self.f4)


    """
    Calculates xn+1 and updates the x value.
    """
    def actualizar_posicion(self):
        self.x += self.h

    """
    Calculates yn+1 and updates the y value (y is equivalent to the velocity).
    """
    def actualizar_velocidad(self):
        self.y += self.h * (25. / 216. * self.f0 + 1408. * self.f2 / 2565. + 2197. *
                            self.f3 / 4104. - self.f4 / 5.)

    """
    Method that calculates the new value of the integration step, according to Fehlberg method.
    If the new h value is smaller than the last one it sets the aux to False, in order to avoid adding
    the same point multiple times, reducing the calculation time. Finally, it updates the value of h.
    """
    def actualizar_paso_integracion(self):
        err = self.h * (1 / 360. * self.f0 - 128 / 4275. * self.f2 - 2197 / 75240. *
                        self.f3 + 1 / 50. * self.f4 + 2 / 55. * self.f5)
        h_nou = 0.9 * self.h * (abs(self.h) * self.error / (abs(err))) ** (1. / 4.)
        if h_nou < self.h:
            self.aux = False      
        self.h = h_nou

    """
    Method that fills lists with the values calculated, in order to graph the results of the calculation.
    More lists could be added if needed.
    """
    def fill_lists(self, lista_x, lista_y, lista_h):
        lista_x.append(self.x)
        lista_y.append(self.y)
        lista_h.append(self.h)




class RungeKutta_2d(RungeKutta):

    """
    Builder of RungeKutta_2d class. We give the needed data in order to solve the problem:
    t: Initial time. It is updated when the calculations are executed.
    x: Initial x position. It is updated when the calculations are executed.
    y: Initial velocity.  It is updated when the calculations are executed.
    xf: Final x position. Given in order to stop the calculation.
    h: Integration step. It may be updated when the calculations are executed if
    step update method (actualizar_paso_integracion) is active.
    error: Desired error.
    This builder works as the father class, but including time to solve problems in two dimensions.
    The methods included are almost equal to the ones in the father class, serve the same purpose and work in the
    same way.
    """
    def __init__(self, t, x, y, tf, h, error, xf=None, *args, **kwargs):
        self.t = t
        self.tf = tf
        self.f = None  #y function
        self.g = None  #x function
        self.aux = True
        """
        With this we apply the heritage from the RungeKutta class
        """
        super(RungeKutta_2d, self).__init__(x, y, xf, h, error, *args, **kwargs)

    """
    In the next methods, calcular_fi, calcular_gi, we implement the calculation of the intermediate functions needed
    for the calculation of xn+1, yn+1 and the update of the integration step. 
    """
    def calcular_f0(self):
        self.f0 = self.f ( self.t,
                           self.x,
                           self.y )
    
    def calcular_f1(self):    
        self.f1 = self.f ( self.t + self.h / 4.,
                           self.x + self.h / 4. * self.g0,
                           self.y + self.h / 4. * self.f0 )
    
    def calcular_f2(self):
        self.f2 = self.f ( self.t + 3. * self.h / 8.,
                           self.x + (3. / 32.) * self.h * self.g0 + (9. / 32.) * self.h * self.g1,
                           self.y + (3. / 32.) * self.h * self.f0 + (9. / 32.) * self.h * self.f1 )
        
    def calcular_f3(self):
        self.f3 = self.f ( self.t + self.h * (12. / 13.),
                           self.x + (1932. / 2197.) * self.h * self.g0 - (7200. / 2197.) * self.h * self.g1 +
                           (7296. / 2197.) * self.h * self.g2,
                           self.y + (1932. / 2197.) * self.h * self.f0 - (7200. / 2197.) * self.h * self.f1 +
                           (7296. / 2197.) * self.h * self.f2 )
        
    def calcular_f4(self):
        self.f4 = self.f ( self.t + self.h,
                           self.x + (439. / 216.) * self.h * self.g0 - 8. * self.h * self.g1 + (3680. / 513.) *
                           self.h * self.g2 - (845. / 4104.) * self.h * self.g3,
                           self.y + (439. / 216.) * self.h * self.f0 - 8. * self.h * self.f1 + (3680. / 513.) *
                           self.h * self.f2 - (845. / 4104.) * self.h * self.f3 )
        
    def calcular_f5(self):
        self.f5 = self.f ( self.t + self.h / 2,
                           self.x - (8. / 27.) * self.h * self.g0 + 2 * self.h * self.g1 - (3544. / 2565.) * self.h *
                           self.g2 + (1859. / 4104.) * self.h * self.g3 - (11. / 40.) * self.h * self.g4,
                           self.y - (8. / 27.) * self.h * self.f0 + 2 * self.h * self.f1 - (3544. / 2565.) * self.h *
                           self.f2 + (1859. / 4104.) * self.h * self.f3 - (11. / 40.) * self.h * self.f4 )
    


    def calcular_g0(self):
        self.g0 = self.g ( self.t, 
                           self.x,
                           self.y )
    
    def calcular_g1(self):
        self.g1 = self.g ( self.t + self.h / 4.,
                           self.x + self.h / 4. * self.g0,
                           self.y + self.h / 4. * self.f0 )
        
    def calcular_g2(self):
        self.g2 = self.g ( self.t + self.h * 3. / 8.,
                           self.x + (3. / 32.) * self.h * self.g0 + (9. / 32.) * self.h * self.g1,
                           self.y + (3. / 32.) * self.h * self.f0 + (9. / 32.) * self.h * self.f1 )
        
    def calcular_g3(self):
        self.g3 = self.g ( self.t + self.h * (12. / 13.),
                           self.x + (1932. / 2197.) * self.h * self.g0 - (7200. / 2197.) * self.h * self.g1 +
                          (7296. / 2197.) * self.h * self.g2, self.y + (1932. / 2197.) * self.h * self.f0 - 
                          (7200. / 2197.) * self.h * self.f1 + (7296. / 2197.) * self.h * self.f2 )
        
    def calcular_g4(self):
        self.g4 = self.g ( self.t + self.h,
                           self.x + (439. / 216.) * self.h * self.g0 - 8. * self.h * self.g1 + (3680. / 513.) *
                          self.h * self.g2 - (845. / 4104.) * self.h * self.g3,
                          self.y + (439. / 216.) * self.h * self.f0 - 8. * self.h * self.f1 + (3680. / 513.) *
                          self.h * self.f2 - (845. / 4104.) * self.h * self.f3 )
        
    def calcular_g5(self):
        self.g5 = self.g( self.t + 0.5 * self.h,
                          self.x - (8. / 27.) * self.h * self.g0 + 2 * self.h * self.g1 - (3544. / 2565.) * self.h *
                          self.g2 + (1859. / 4104.) * self.h * self.g3 - (11. / 40.) * self.h * self.g4,
                          self.y - (8. / 27.) * self.h * self.f0 + 2 * self.h * self.f1 - (3544. / 2565.) * self.h *
                          self.f2 + (1859. / 4104.) * self.h * self.f3 - (11. / 40.) * self.h * self.f4 )


    
    def actualizar_tiempo(self):
        self.t = self.t + self.h
        
    def actualizar_posicion(self):
        self.x = self.x + self.h * (25. / 216. * self.g0 + 1408. * self.g2 / 2565. + 2197. * self.g3 / 4104. - self.g4 / 5.)
    
    def actualizar_velocidad(self):
        self.y = self.y + self.h * (25. / 216. * self.f0 + 1408. * self.f2 / 2565. + 2197. * self.f3 / 4104. - self.f4 / 5.)


    """
    In the method that calculates the new value of the integration step we calculate one for each function (h_nou_y and 
    h_nou_x) and then we take the smallest of the two, in order to diminish errors.
    """
    def actualizar_paso_integracion(self):
        Err1 = self.h * (1 / 360. * self.g0 - 128 / 4275. * self.g2 - 2197 / 75240. * self.g3 + 1 / 50. * self.g4 + 2 / 55. * self.g5)
        Err2 = self.h * (1 / 360. * self.f0 - 128 / 4275. * self.f2 - 2197 / 75240. * self.f3 + 1 / 50. * self.f4 + 2 / 55. * self.f5)
        h_nou_y = 0.9 * self.h * (self.h * self.error / (abs (Err1))) ** (1. / 4.)
        h_nou_x = 0.9 * self.h * (self.h * self.error / (abs (Err2))) ** (1. / 4.)
        h_nou = min (h_nou_x, h_nou_y)

        if h_nou < self.h:
            self.r = False
            self.h = h_nou

    def fill_lists(self, lista_t, lista_x, lista_y, lista_h):
        lista_t.append(self.t)
        lista_x.append(self.x)
        lista_y.append(self.y)
        lista_h.append(self.h)