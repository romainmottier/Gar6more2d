import pandas as pd
import csv
import math 
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
# import scipy as sc
from scipy import signal


####################################################################################################
####################################################################################################
####################################################################################################

plt.close()
# Obtenez le chemin absolu du fichier en cours
current_file_path = os.path.abspath(__file__)
# Obtenez le répertoire du fichier en cours
current_directory = os.path.dirname(current_file_path)
# Changez le répertoire de travail actuel pour le répertoire du fichier en cours
os.chdir(current_directory)
 
# Sensors
k1 = np.loadtxt('Interface_k1.csv', delimiter = ',', skiprows = 1, dtype = float)
k2 = np.loadtxt('Interface_k2.csv', delimiter = ',', skiprows = 1, dtype = float)
k3 = np.loadtxt('Interface_k3.csv', delimiter = ',', skiprows = 1, dtype = float)
k4 = np.loadtxt('Interface_k4.csv', delimiter = ',', skiprows = 1, dtype = float)

# Variables 
dt = 0.1*2**(-9)
Tf = 0.25
t  = Tf*k1[:,0]*dt
delay = 0.15
# t  = t + delay 
k1 = k1[:,1]
k2 = k2[:,1]
k3 = k3[:,1]
k4 = k4[:,1]

k1 = k1*3500
plt.plot(t, k1, 
            linestyle = "-", 
            linewidth = 2,
            label     = "$\mathbf{k = 1}$",
            color     = "darkgreen")

k2 = k2*3500
plt.plot(t, k2, 
            linestyle = "-", 
            linewidth = 2,
            label     = "$\mathbf{k = 2}$",
            color     = "orange")

k3 = k3*3500
plt.plot(t, k3, 
            linestyle = "-", 
            linewidth = 2,
            label     = "$\mathbf{k = 3}$",
            color     = "darkred")

k4 = k4*3500
plt.plot(t, k4, 
            linestyle = "-", 
            linewidth = 2,
            label     = "$\mathbf{k = 4}$",
            color     = "navy")

# Gar6more
Ux = pd.read_csv("Ux.dat", sep='\s+', usecols=[0, 1], names=["T", "Ux"])
Uy = pd.read_csv("Uy.dat", sep='\s+', usecols=[0, 1], names=["T", "Uy"])

def lim_pente_1(y):
     y = np.array(y)
     n = len(y)
     ener=1
     for dx in range(1,100):
         y0=y.copy()
         for i in range(dx,n-dx):
             ymin=min(y0[i-dx],y0[i+dx])
             ymax=max(y0[i-dx],y0[i+dx])
             y[i]=max(min(y0[i],ymax),ymin)
     return y
def lim_pente_2(y):
    y = np.array(y)
    nx =len(y)
    y0 = np.zeros(nx)
    ener=1
    for dx in range(1,10):
        y0[0:nx]=y[0:nx]   
        for i in range(dx,nx-dx):
            ymin=min(y0[i-dx],y0[i+dx])
            ymax=max(y0[i-dx],y0[i+dx])
            y[i]=max(min(y0[i],ymax),ymin)
        ener0=ener            
        ener= np.sqrt((y-y0).dot(     y-y0))
        if(dx==1):
            e0=ener
        ener/=e0
        if(dx>1 and ener>ener0):
            y=y0
            break
    return y
Ux["Ux"] = lim_pente_1(Ux["Ux"])
Ux["Ux"] = signal.savgol_filter(np.array(Ux["Ux"]),5,2)


Ux["T"] = Ux["T"] - delay
plt.plot(Ux["T"].to_numpy(), -Ux["Ux"].to_numpy(), 
            linestyle="-", 
            linewidth=2,
            label="$\mathbf{Gar6more}$",
            color="k")

#################### Legends
plt.legend(fontsize="12.5", framealpha=0) # Legends
plt.xlabel("Time")                        # Abscissas title
plt.ylabel("$v_x$")   
plt.grid()
plt.show()


plt.savefig('academic_pulse_vy.png', transparent=True)
