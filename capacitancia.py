import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------
# Parámetros
# ----------------------------------------------------
W=L=85e-6

# ----------------------------------------------------
# Función para leer un archivo tipo SPICE
# ----------------------------------------------------
def leer_data(filename):
    data = np.loadtxt(filename, skiprows=1)  # saltea el header
    t  = data[:,0]
    vg = data[:,1]
    ig = data[:,2]
    return t, vg, ig

# ----------------------------------------------------
# Leer ambos archivos
# ----------------------------------------------------
t1, vg1, ig1 = leer_data("data_1.txt")  # Teórico
t2, vg2, ig2 = leer_data("data_2.txt")  # Real

# ----------------------------------------------------
# Calcular capacitancia
# ----------------------------------------------------
C1 = (1/(W * L)) * ig1 #en función de Vg1 
C2 = (1/(W * L)) * ig2 #en función de vg2

# ----------------------------------------------------
# Gráfico
# ----------------------------------------------------
plt.figure(figsize=(8,5))
plt.plot(vg1, C1, label="Teórico")
plt.plot(vg2, C2, label="Real", linestyle='--')

plt.xlabel("Vg (V)")
plt.ylabel("Capacitancia (F)")
plt.title("Capacitancia vs Vg")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
