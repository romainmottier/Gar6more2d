# import pandas as pd
# import csv
# import math 
# import os
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MaxNLocator

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################

# plt.close()
# # Obtenez le chemin absolu du fichier en cours
# current_file_path = os.path.abspath(__file__)
# # Obtenez le répertoire du fichier en cours
# current_directory = os.path.dirname(current_file_path)
# # Changez le répertoire de travail actuel pour le répertoire du fichier en cours
# os.chdir(current_directory)
 
# # Sensors
# k1 = np.loadtxt('Acoustic_k1.csv', delimiter = ',', skiprows = 1, dtype = float)
# k2 = np.loadtxt('Acoustic_k2.csv', delimiter = ',', skiprows = 1, dtype = float)
# k3 = np.loadtxt('Acoustic_k3.csv', delimiter = ',', skiprows = 1, dtype = float)
# k4 = np.loadtxt('Acoustic_k4.csv', delimiter = ',', skiprows = 1, dtype = float)

# # Variables 
# dt = 0.1*2**(-9)
# Tf = 0.25
# t  = Tf*k1[:,0]*dt
# delay = 0.15
# # t  = t + delay 
# k1 = k1[:,1]
# k2 = k2[:,1]
# k3 = k3[:,1]
# k4 = k4[:,1]

# k1 = k1*3500
# plt.plot(t, k1, 
#             linestyle = "-", 
#             linewidth = 2,
#             label     = "$\mathbf{k = 1}$",
#             color     = "darkgreen")

# k2 = k2*3500
# plt.plot(t, k2, 
#             linestyle = "-", 
#             linewidth = 2,
#             label     = "$\mathbf{k = 2}$",
#             color     = "orange")

# k3 = k3*3500
# plt.plot(t, k3, 
#             linestyle = "-", 
#             linewidth = 2,
#             label     = "$\mathbf{k = 3}$",
#             color     = "darkred")

# k4 = k4*3500
# plt.plot(t, k4, 
#             linestyle = "-", 
#             linewidth = 2,
#             label     = "$\mathbf{k = 4}$",
#             color     = "navy")

# # Gar6more
# P = pd.read_csv("P.dat", delim_whitespace=True, usecols=[0, 1], names=["T", "P"])

# P["T"] = P["T"] - delay
# plt.plot(P["T"].to_numpy(), P["P"].to_numpy(), 
#             linestyle="-", 
#             linewidth=2,
#             label="$\mathbf{Gar6more}$",
#             color="k")

# #################### Legends
# plt.gca().spines['top'].set_linewidth(1.5)    # Bordure supérieure
# plt.gca().spines['bottom'].set_linewidth(1.5) # Bordure inférieure
# plt.gca().spines['left'].set_linewidth(1.5)   # Bordure gauche
# plt.gca().spines['right'].set_linewidth(1.5)  # Bordure droite

# plt.legend(fontsize="20", framealpha=1, loc='upper left')                        # Legends
# plt.xlabel("Time", fontsize=20, fontweight='bold')                     # Abscissas title
# plt.ylabel("Acoustic pressure", fontsize=20, fontweight='bold')   
# plt.grid(True, linewidth=2)
# plt.show()

# plt.savefig('academic_pulse_pressure.png', transparent=True)


# import pandas as pd
# import csv
# import math 
# import os
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MaxNLocator

# plt.close()
# # Obtenez le chemin absolu du fichier en cours
# current_file_path = os.path.abspath(__file__)
# # Obtenez le répertoire du fichier en cours
# current_directory = os.path.dirname(current_file_path)
# # Changez le répertoire de travail actuel pour le répertoire du fichier en cours
# os.chdir(current_directory)
 
# # Sensors
# k1 = np.loadtxt('Acoustic_k1.csv', delimiter = ',', skiprows = 1, dtype = float)
# k2 = np.loadtxt('Acoustic_k2.csv', delimiter = ',', skiprows = 1, dtype = float)
# k3 = np.loadtxt('Acoustic_k3.csv', delimiter = ',', skiprows = 1, dtype = float)
# k4 = np.loadtxt('Acoustic_k4.csv', delimiter = ',', skiprows = 1, dtype = float)

# # Variables 
# dt = 0.1*2**(-9)
# Tf = 0.25
# t  = Tf*k1[:,0]*dt
# delay = 0.15
# # t  = t + delay 
# k1 = k1[:,1]
# k2 = k2[:,1]
# k3 = k3[:,1]
# k4 = k4[:,1]

# k1 = k1*3500
# k2 = k2*3500
# k3 = k3*3500
# k4 = k4*3500

# # Gar6more
# P = pd.read_csv("P.dat", delim_whitespace=True, usecols=[0, 1], names=["T", "P"])
# P["T"] = P["T"] - delay

# # Création de la figure et de la grille
# fig, ax = plt.subplots(figsize=(20, 15))

# # Graphique principal
# ax.plot(t, k1, linestyle="-", linewidth=2, label="$\mathbf{k = 1}$", color="darkgreen")
# ax.plot(t, k2, linestyle="-", linewidth=2, label="$\mathbf{k = 2}$", color="orange")
# ax.plot(t, k3, linestyle="-", linewidth=2, label="$\mathbf{k = 3}$", color="darkred")
# ax.plot(t, k4, linestyle="-", linewidth=2, label="$\mathbf{k = 4}$", color="navy")
# ax.plot(P["T"].to_numpy(), P["P"].to_numpy(), linestyle="-", linewidth=2, label="$\mathbf{Gar6more}$", color="k")

# # Personnalisation des bordures
# ax.spines['top'].set_linewidth(1.5)
# ax.spines['bottom'].set_linewidth(1.5)
# ax.spines['left'].set_linewidth(1.5)
# ax.spines['right'].set_linewidth(1.5)

# # Légende
# legend = ax.legend(fontsize="20", framealpha=1, loc='upper left')

# # Annotation ou schéma
# # Vous pouvez ajuster les coordonnées et le texte selon vos besoins
# ax.annotate('Schéma ici', xy=(0.25, 0.85), xycoords='axes fraction', fontsize=14,
#             bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))

# ax.set_xlabel("Time", fontsize=20, fontweight='bold')
# ax.set_ylabel("Acoustic pressure", fontsize=20, fontweight='bold')
# ax.grid(True, linewidth=2)

# # Affichage du graphique
# plt.show()

# # Sauvegarde du graphique
# fig.savefig('academic_pulse_pressure.png', transparent=True)


# import pandas as pd
# import csv
# import math 
# import os
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MaxNLocator
# from matplotlib.patches import Rectangle

# plt.close()
# # Obtenez le chemin absolu du fichier en cours
# current_file_path = os.path.abspath(__file__)
# # Obtenez le répertoire du fichier en cours
# current_directory = os.path.dirname(current_file_path)
# # Changez le répertoire de travail actuel pour le répertoire du fichier en cours
# os.chdir(current_directory)
 
# # Sensors
# k1 = np.loadtxt('Acoustic_k1.csv', delimiter = ',', skiprows = 1, dtype = float)
# k2 = np.loadtxt('Acoustic_k2.csv', delimiter = ',', skiprows = 1, dtype = float)
# k3 = np.loadtxt('Acoustic_k3.csv', delimiter = ',', skiprows = 1, dtype = float)
# k4 = np.loadtxt('Acoustic_k4.csv', delimiter = ',', skiprows = 1, dtype = float)

# # Variables 
# dt = 0.1*2**(-9)
# Tf = 0.25
# t  = Tf*k1[:,0]*dt
# delay = 0.15
# # t  = t + delay 
# k1 = k1[:,1]
# k2 = k2[:,1]
# k3 = k3[:,1]
# k4 = k4[:,1]

# k1 = k1*3500
# k2 = k2*3500
# k3 = k3*3500
# k4 = k4*3500

# # Gar6more
# P = pd.read_csv("P.dat", delim_whitespace=True, usecols=[0, 1], names=["T", "P"])
# P["T"] = P["T"] - delay

# # Création de la figure et de la grille
# fig, ax = plt.subplots(figsize=(10, 6))

# # Graphique principal
# ax.plot(t, k1, linestyle="-", linewidth=2, label="$\mathbf{k = 1}$", color="darkgreen")
# ax.plot(t, k2, linestyle="-", linewidth=2, label="$\mathbf{k = 2}$", color="orange")
# ax.plot(t, k3, linestyle="-", linewidth=2, label="$\mathbf{k = 3}$", color="darkred")
# ax.plot(t, k4, linestyle="-", linewidth=2, label="$\mathbf{k = 4}$", color="navy")
# ax.plot(P["T"].to_numpy(), P["P"].to_numpy(), linestyle="-", linewidth=2, label="$\mathbf{Gar6more}$", color="k")

# # Personnalisation des bordures
# ax.spines['top'].set_linewidth(1.5)
# ax.spines['bottom'].set_linewidth(1.5)
# ax.spines['left'].set_linewidth(1.5)
# ax.spines['right'].set_linewidth(1.5)

# # Légende
# legend = ax.legend(fontsize="20", framealpha=1, loc='upper left')

# # Dessiner un carré à côté de la légende
# square_x = 0.7
# square_y = 0.85
# square_size = 0.1

# # Conversion des coordonnées de fraction de la boîte englobante de l'axe
# bbox = ax.get_position()
# fig_width, fig_height = fig.get_size_inches()

# square_x_fig = bbox.x0 + square_x * (bbox.x1 - bbox.x0)
# square_y_fig = bbox.y0 + square_y * (bbox.y1 - bbox.y0)

# rect = Rectangle((square_x_fig * fig_width, square_y_fig * fig_height), 
#                  square_size * fig_width, 
#                  square_size * fig_height,
#                  linewidth=1, edgecolor='black', facecolor='none', transform=fig.transFigure)

# fig.patches.append(rect)

# ax.set_xlabel("Time", fontsize=20, fontweight='bold')
# ax.set_ylabel("Acoustic pressure", fontsize=20, fontweight='bold')
# ax.grid(True, linewidth=2)

# # Affichage du graphique
# plt.show()

# # Sauvegarde du graphique
# fig.savefig('academic_pulse_pressure.png', transparent=True)


import pandas as pd
import csv
import math 
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle

plt.close()
# Obtenez le chemin absolu du fichier en cours
current_file_path = os.path.abspath(__file__)
# Obtenez le répertoire du fichier en cours
current_directory = os.path.dirname(current_file_path)
# Changez le répertoire de travail actuel pour le répertoire du fichier en cours
os.chdir(current_directory)
 
# Sensors
k1 = np.loadtxt('Acoustic_k1.csv', delimiter = ',', skiprows = 1, dtype = float)
k2 = np.loadtxt('Acoustic_k2.csv', delimiter = ',', skiprows = 1, dtype = float)
k3 = np.loadtxt('Acoustic_k3.csv', delimiter = ',', skiprows = 1, dtype = float)
k4 = np.loadtxt('Acoustic_k4.csv', delimiter = ',', skiprows = 1, dtype = float)

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
k2 = k2*3500
k3 = k3*3500
k4 = k4*3500

# Gar6more
P = pd.read_csv("P.dat", delim_whitespace=True, usecols=[0, 1], names=["T", "P"])
P["T"] = P["T"] - delay

# Création de la figure et de la grille
fig, ax = plt.subplots(figsize=(10, 6))

# Graphique principal
ax.plot(t, k1, linestyle="-", linewidth=2, label="$\mathbf{k = 1}$", color="darkgreen")
ax.plot(t, k2, linestyle="-", linewidth=2, label="$\mathbf{k = 2}$", color="orange")
ax.plot(t, k3, linestyle="-", linewidth=2, label="$\mathbf{k = 3}$", color="darkred")
ax.plot(t, k4, linestyle="-", linewidth=2, label="$\mathbf{k = 4}$", color="navy")
ax.plot(P["T"].to_numpy(), P["P"].to_numpy(), linestyle="-", linewidth=2, label="$\mathbf{Gar6more}$", color="k")

# Personnalisation des bordures
ax.spines['top'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
ax.spines['right'].set_linewidth(1.5)

# Légende
legend = ax.legend(fontsize="20", framealpha=1, loc='upper left')

# Dessiner un carré à côté de la légende
square_x = 0.7
square_y = 0.85
square_size = 0.1

# Conversion des coordonnées de fraction de la boîte englobante de l'axe
bbox = ax.get_position()
fig_width, fig_height = fig.get_size_inches()

square_x_fig = bbox.x0 + square_x * (bbox.x1 - bbox.x0)
square_y_fig = bbox.y0 + square_y * (bbox.y1 - bbox.y0)

rect = Rectangle((square_x_fig * fig_width, square_y_fig * fig_height), 
                 square_size * fig_width, 
                 square_size * fig_height,
                 linewidth=1, edgecolor='black', facecolor='none', transform=fig.transFigure)

fig.patches.append(rect)

ax.set_xlabel("Time", fontsize=20, fontweight='bold')
ax.set_ylabel("Acoustic pressure", fontsize=20, fontweight='bold')
ax.grid(True, linewidth=2)

# Affichage du graphique
plt.show()

# Sauvegarde du graphique
fig.savefig('academic_pulse_pressure.png', transparent=True)
