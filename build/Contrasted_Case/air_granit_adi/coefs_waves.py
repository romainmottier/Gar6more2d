
import numpy as np

def lame_coefficients_solid(Vp, Vs, rho):
    mu = rho * Vs**2
    lambda_ = rho*Vp**2 -2*mu
    return lambda_, mu

def lame_coefficients_fluid(Vp2, rho):
    mu = rho * Vp2**2
    return mu

# Exemple d'utilisation :
Vp_s  = 17.5  # Vitesse des ondes P en m/s
Vs    = 9.0          # Vitesse des ondes S en m/s
rho_s = 2200         # Masse volumique en kg/m^3
Vp_f  = 1.0
rho_f = 1.0

lambda_, mu = lame_coefficients_solid(Vp_s, Vs, rho_s)
print("\n")
print("SOLIDE:")
print("mu =", mu)
print("lambda =", lambda_)
print("rho =", rho_s)
print("\n")
mu = lame_coefficients_fluid(Vp_f, rho_f)
print("FLUIDE:")
print("mu =", mu)
print("rho =", rho_f)
print("\n")

