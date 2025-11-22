import numpy as np
import matplotlib.pyplot as plt

#Parametros
e_o = 8.854e-12 #F/m
epsilon_si = 11.7*e_o
epsilon_ox = 3.9*e_o


t_ox = 16e-9 #m 

q = 1.602e-19       # Carga elemental (C)
k = 1.381e-23       # Constante de Boltzmann (J/K)
T = 300             # Temperatura (K)

V_th = (k*T)/q
