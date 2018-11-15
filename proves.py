#!/usr/bin/env python
# coding: utf-8

# In[3]:



import numpy as np
import matplotlib.pyplot as plt

from collections import namedtuple


# In[4]:


PosicionTiempo = namedtuple('PosicionTiempo', ('posicion', 'tiempo')) 


# In[5]:


def vert_shot(x_0, v_0, dt, g = -9.81):
    x = x_0
    v = v_0
    i = t = 0
    pos = []
    temps = []
    while x >= 0:
        i = i+1
        pos.append(x)
        temps.append(t)
        x = x + v * dt
        v = v + g * dt
        t = i * dt
    pos.append(x)
    temps.append(t)   
    
    return PosicionTiempo(pos, temps)


# In[6]:


# Tiro vertical forward in time

def vert_shot_ft(x_0, v_0, dt, g = -9.81):
    x = x_0
    v = v_0
    i = t = 0
    pos = []
    temps = []
    while x >= 0:
        i = i+1
        pos.append(x)
        temps.append(t)
        v = v + g * dt
        x = x + v * dt
        t = i * dt
    pos.append(x)
    temps.append(t)   
    
    return PosicionTiempo(pos, temps)


# In[7]:


# Tiro vertical centered in time
def vert_shot_ct(x_0, v_0, dt, g = -9.81):
    x = x_0
    v = v_0
    i = t = 0
    pos = []
    temps = []
    while x >= 0:
        i = i+1
        pos.append(x)
        temps.append(t)
        x = x + (v + 0.5 * g * dt) * dt
        v = v + g * dt
        t = i * dt
    pos.append(x)
    temps.append(t) 
    return PosicionTiempo(pos, temps)


# In[8]:


def vert_shot_analytic(x_0, v_0, dt, g = -9.81):
    x = x_0
    i = t = 0
    pos = []
    temps = []
    pos.append(x)
    temps.append(t)
    while x >= 0:
        i = i+1
        t = i * dt
        x = v_0 * t + 0.5 * g * t**2
        pos.append(x)
        temps.append(t)
    
    return PosicionTiempo(pos, temps)


# In[9]:


x_0 = 0
v_0 = 25
dt = 0.1
dt_2 = 0.5
dt_3 = 0.01


# In[10]:


pt1 = vert_shot(x_0, v_0, dt)


# In[11]:


pt2 = vert_shot_ft(x_0, v_0, dt)


# In[12]:


pt3 = vert_shot_ct(x_0, v_0, dt)


# In[13]:


pta = vert_shot_analytic(x_0, v_0, dt)


# In[17]:


pt1_2 = vert_shot(x_0, v_0, dt_2)


# In[18]:


pt2_2 = vert_shot_ft(x_0, v_0, dt_2)


# In[19]:


pt3_2 = vert_shot_ct(x_0, v_0, dt_2)


# In[20]:


pt1_3 = vert_shot(x_0, v_0, dt_3)


# In[21]:


pt2_3 = vert_shot_ft(x_0, v_0, dt_3)


# In[22]:


pt3_3 = vert_shot_ct(x_0, v_0, dt_3)


# In[23]:


plt.title("Vertical Shot, dt = 0.1")
plt.xlabel("Time")
plt.ylabel("Position")
plt.plot(pt1.tiempo, pt1.posicion, label="Numerical Calc. BT")
plt.plot(pt2.tiempo, pt2.posicion, label="Numerical Calc. FT")
plt.plot(pt3.tiempo, pt3.posicion, label="Numerical Calc. CT")
plt.plot(pta.tiempo, pta.posicion, label = "Analytical Solution")
plt.legend()
plt.show()

# In[24]:


plt.title("Vertical Shot, dt = 0.5")
plt.xlabel("Time")
plt.ylabel("Position")
plt.plot(pt1_2.tiempo, pt1_2.posicion, label="Numerical Calc. BT")
plt.plot(pt2_2.tiempo, pt2_2.posicion, label="Numerical Calc. FT")
plt.plot(pt3_2.tiempo, pt3_2.posicion, label="Numerical Calc. CT")
plt.plot(pta.tiempo, pta.posicion, label = "Analytical Solution")
plt.legend()
plt.show()

# In[25]:


plt.title("Vertical Shot, dt = 0.01")
plt.xlabel("Time")
plt.ylabel("Position")
plt.plot(pt1_3.tiempo, pt1_3.posicion, label="Numerical Calc. BT")
plt.plot(pt2_3.tiempo, pt2_3.posicion, label="Numerical Calc. FT")
plt.plot(pt3_3.tiempo, pt3_3.posicion, label="Numerical Calc. CT")
plt.plot(pta.tiempo, pta.posicion, label = "Analytical Solution")
plt.legend()
plt.show()

# In[ ]:




