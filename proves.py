#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import matplotlib.pyplot as plt


# In[4]:


# Cálculo del número pi

def calc_pi(n):
    for i in range (3, n):
        x = i * np.sin(np.pi/i)
    return x, np.pi-x

print calc_pi(100000)


# In[5]:


# Tiro vertical
def vert_shot(x_0, v_0, dt, g = -9.81):
    x = x_0
    v = v_0
    i = 0
    pos = []
    temps = []
    while x >= 0:
        i = i+1
        x_p = x + v * dt
        v_p = v + g * dt
        v = v_p
        x = x_p
        t = i * dt
        pos.append(x)
        temps.append(t)
        print (x, v)
    
    plt.plot(temps, pos)
    plt.show()
        

vert_shot(0, 25, 0.01)
        
        


# In[ ]:




