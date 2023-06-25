import numpy as np
import matplotlib.pyplot as plt
from numba import njit


##Parameter definitions (SI)
alphaConst = 30
theta = 70
h = 0.4 # m
l = 0.7 # m
a = 0.55 # m
b = 0.15 # m
c = l/2 # m
d = h/2 # m
e = 0.03 # m
f = 2000 # N
g = 0.18 # m
i = 0.25 # m

@njit
def F(t): # N (SI)
  return f + (f/2)*np.sin(t/2) + (f/20)*np.sin(20*t) 

@njit
def alpha(w, amortecimento): # Atencao! O vetor w deve estar em radianos!
  return (2*amortecimento*w[0]*w[5])/(w[0] + w[5]) 

@njit
def beta(w, amortecimento): # Atencao! O vetor w deve estar em radianos!
  return (2*amortecimento)/(w[0] + w[5]) 

@njit
def modulo_u(ux, uy): # Atencao! O vetor w deve estar em radianos!
  return np.sqrt(ux**2 + uy**2) 


def main():
    # a1
    delta_x = 0.5
    delta_x = 1
    delta_x = 2

    # a2
    amortecimento = 0.06
    t_final = 10 # s (SI)


if __name__ == "__main__":
    main()
