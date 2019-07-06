# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 20:37:35 2017

@author: Michael
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import random

data = np.loadtxt("csim__max_target-displace-sphere_10.csv", delimiter = ",")
new_data = data.flatten(order = 'C')
q = np.arange(0,1000,1)
const_array = np.array([])
mean_array = np.zeros(1000)

def exp_decay(t,const):
    this = 82492542718.918503*np.exp(-const*t)
    return this
    
def time_decay(const,nt):
    this = np.log(nt/82492542718.918503)/-const
    return this
    
sorted_data = np.sort(new_data)
    
    
condition = (new_data < 1.5E9)
    
smaller_data = np.extract(condition, new_data)
print(len(smaller_data))
print(smaller_data)
print((len(smaller_data)/len(new_data))*100)
print(min(smaller_data))
print(len(new_data))
    
min_array = np.zeros_like(new_data)
n = np.arange(1,len(new_data),1)
min_point = new_data[0]
min_array[0] = new_data[0]
    
drop_points = np.array([new_data[0]])
drop_axis = np.array([0])
for i in n:
    if new_data[i] < min_point:
        min_point = new_data[i]
        drop_points = np.append(drop_points,new_data[i])
        drop_axis = np.append(drop_axis,i)
        min_array[i] = min_point
    

plt.figure()
plt.hist(new_data, 100)
plt.xlabel("Closest distance to earth [m]")
plt.ylabel("Number of Ejecta")

"""


t = time_decay(0.317,(6000E3))
#x = exp_decay(t,0.317)

x = np.arange(0,200000,1)
y = exp_decay(x,0.317)
plt.plot(x,y)

"""

"""
data = np.loadtxt("csim__max_target-displace-sphere_2000.csv", delimiter = ",")
new_data = data.flatten(order = 'A')
q = np.arange(0,1000,1)
const_array = np.array([])
mean_array = np.zeros(1000)

for j in q:
    random.shuffle(new_data)
    
    sorted_data = np.sort(new_data)
    
    
    condition = (new_data < 1.5E9)
    
    smaller_data = np.extract(condition, new_data)
    #print(len(smaller_data))
    #print(smaller_data)
    #print((len(smaller_data)/len(new_data))*100)
   # print(min(smaller_data))
   # print(len(new_data))
    
    min_array = np.zeros_like(new_data)
    n = np.arange(1,len(new_data),1)
    min_point = new_data[0]
    min_array[0] = new_data[0]
    
    drop_points = np.array([new_data[0]])
    drop_axis = np.array([0])
    for i in n:
        if new_data[i] < min_point:
            min_point = new_data[i]
            drop_points = np.append(drop_points,new_data[i])
            drop_axis = np.append(drop_axis,i)
        min_array[i] = min_point
    
    
    #plt.figure()     
   # plt.plot(drop_axis,drop_points, "x") 
         
    
    popt, pcov = so.curve_fit(exp_decay,drop_axis, drop_points, p0 = (0.001))
    x = np.arange(0,200000,1)
    y = exp_decay(x,popt[0])
    #plt.plot(x,y)
    
    #plt.figure()
    #plt.plot(sorted_data)
    
    const_array = np.append(const_array, popt[0])
    mean_array[j] = np.mean(const_array)
    #print(popt[0])
    print(j)

print(np.mean(const_array))

plt.figure()
plt.plot(mean_array)

"""
