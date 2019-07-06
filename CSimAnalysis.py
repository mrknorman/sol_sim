# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 17:06:46 2017

@author: Michael Norman

"""

""" ~~~~~~~~~~ Importing ~~~~~~~~~~ """  

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D 
import numpy as np

""" ~~~~~~~~~~ Functions ~~~~~~~~~~ """  

def formatGraph(g_title,x_title,y_title):
    
    """ Ease of use function """
    
    plt.title(g_title) 
    plt.xlabel(x_title) 
    plt.legend(loc='best')
    plt.ylabel(y_title) 
    plt.grid()

no_planets = 6
no_ejecta = 36

no_stuff = (no_ejecta + no_planets)*3

data = np.loadtxt("csim_posit_12.csv", delimiter = ",", unpack = True, usecols = (range(no_stuff)))

object_name = ["Sun","Mercury","Venus","Earth", "Mars","Jupiter"]

x = 0

plt.figure()
for i in object_name:
    plt.plot(data[x],data[x+1], label = i)
    
    x += 3
   
n = np.arange(no_planets *3,(no_planets+no_ejecta)*3,3)

for i in n:
    plt.plot(data[i],data[i+1])
    
formatGraph("Solar System Sim","X[m]","Y[m]")

"""

data2 = np.loadtxt("csim_origin_displace_1.csv", delimiter = ",", unpack = True, usecols = (range(no_ejecta)))

n = np.arange(0,no_ejecta,1)
aphellion = np.zeros(no_ejecta)

for i in n:

    aphellion[i] = max(data2[i])

plt.figure()
plt.plot(aphellion)
"""

   
""" ~~~~~~~~~~ 3d plot  ~~~~~~~~~~ """  


"""
fig1 = plt.figure()
ax = Axes3D(fig1)
max_range = np.array([])

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

n = np.arange(0,7*3,3)

for i in n:

    theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
    z = np.array(data[n + 2])  
    r = z**2 + 1
    x = np.array(data[n + 0])
    y = np.array(data[n + 1])
    ax.plot(x, y, zs = z, label='parametric curve')
    ax.legend()

max_range = np.append(max_range, x.max())
max_range = np.append(max_range, y.max())
max_range = np.append(max_range, z.max())
                 
MaxRange = max_range.max()
fig.set_xlim(- MaxRange, + MaxRange)
fig.set_ylim(- MaxRange, + MaxRange)
fig.set_zlim(- MaxRange, + MaxRange)
        

    
    
plt.show()

for i in n:
    X = np.array(data[n + 0])
    Y = np.array(data[n + 1])
    Z = np.array(data[n + 2])  
    ax.plot(X,Y,Z)

    
"""