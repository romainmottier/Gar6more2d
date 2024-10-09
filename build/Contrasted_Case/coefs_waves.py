
import numpy as np

def lame_coefficients_solid(Vp, Vs, rho):
    mu = rho * Vs**2
    lambda_ = rho*Vp**2 -2*mu
    return lambda_, mu

def lame_coefficients_fluid(Vp2, rho):
    mu = rho * Vp2**2
    return mu

# Exemple d'utilisation :
Vp_s  = 6000.0  # Vitesse des ondes P en m/s
Vs    = 3000.0          # Vitesse des ondes S en m/s
rho_s = 2690.0         # Masse volumique en kg/m^3
Vp_f  = 1500.0
rho_f = 1025.0

lambda_, mu = lame_coefficients_solid(Vp_s, Vs, rho_s)
print("\n")
print("SOLIDE:")
print("lambda =", lambda_)
print("mu =", mu)
print("\n")
mu = lame_coefficients_fluid(Vp_f, rho_f)
print("FLUIDE:")
print("mu =", mu)
print("\n")

