# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 16:27:13 2017

@author: Michael
"""

import numpy as np
import matplotlib.pyplot as plt


no_ejecta = 36
Blank = np.loadtxt("csim__max_target-displace0-50_1.csv", delimiter = ",", unpack = True, skiprows = 0, usecols = (range(no_ejecta)))

xlabel = np.linspace(0,360,20)
ylabel = np.linspace(0,360,36)

plt.figure()
plt.contourf(xlabel,ylabel,Blank)
col = plt.colorbar()
col.set_label("Distance to Earth[m]", fontsize=20)
#plt.contour(xlabel,ylabel,Blank, linewidths = 4)

plt.xlabel("Angle around star [Deg]", fontsize=20)
plt.ylabel("Angle around planet [Deg]", fontsize=20)


