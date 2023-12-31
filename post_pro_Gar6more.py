import pandas as pd
import matplotlib.pyplot as plt

UX = pd.read_csv("Ux.dat", delim_whitespace=True, usecols=[0, 1], names=["T", "Ux"])
UY = pd.read_csv("Uy.dat", delim_whitespace=True, usecols=[0, 1], names=["T", "Uy"])
P = pd.read_csv("P.dat", delim_whitespace=True, usecols=[0, 1], names=["T", "P"])

fig = plt.figure(figsize=(10, 5))
gs = fig.add_gridspec(1, 3, wspace=0.25)
axs = gs.subplots()

axs[0].plot(UX["T"], UX["Ux"])
axs[1].plot(UY["T"], UY["Uy"])
axs[2].plot(P["T"], P["P"])

plt.show()
