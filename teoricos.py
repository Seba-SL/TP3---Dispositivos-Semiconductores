# -------------------------------------------------------------
# Cálculo de curvas teóricas para el Grupo 6 de la tabla MOS
# Qs'(psi_s), psi_s(VG)
# -------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# Constantes físicas
q = 1.602176634e-19       # C
eps0 = 8.8541878128e-12   # F/m
k_B = 1.380649e-23        # J/K
T = 300.0                 # K
Vth = k_B * T / q         # ~25.85 mV

# Parámetros del Grupo 6
t_ox = 16e-9              # 16 nm
N_A_cm3 = 2e16            # cm^-3 (dopado p)
N_A = N_A_cm3 * 1e6       # m^-3
eps_ox = eps0 * 3.9       # SiO2
eps_si = eps0 * 11.7      # Si

Cox = eps_ox / t_ox       # F/m^2

# Suposiciones típicas (modificar si las conocés)
V_FB = -1.0               # Voltaje de flat-band ASSUMIDO
Q_ox = 0.0                # Sin carga en óxido

# Rango de ψs
psi = np.linspace(1e-6, 1.2, 1000)

# Carga en la superficie en deplexión:
Qs = -np.sqrt(2 * q * eps_si * N_A * psi)

# Ecuación MOS:
# VG = V_FB + psi - Qs / Cox
V_G = V_FB + psi - Qs / Cox

# Barrido de VG
V_G_sweep = np.linspace(-2.0, 1.0, 1000)
psi_from_VG = np.interp(V_G_sweep, V_G, psi, left=np.nan, right=np.nan)

# --------- GRÁFICOS ----------
plt.figure(figsize=(6,4))
plt.plot(psi, Qs, alpha = 0.9, linewidth=3)
plt.xlabel(r'$\psi_s$ (V)')
plt.ylabel(r"$Q'_s$ (C/m$^2$)")
plt.title("Grupo 6: $Q'_s$ vs $\\psi_s$")
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,4))
plt.plot(V_G, psi,alpha = 0.9, linewidth=3)
plt.xlabel(r'$V_G$ (V)')
plt.ylabel(r'$\psi_s$ (V)')
plt.title("Grupo 6: $\\psi_s$ vs $V_G$")
plt.grid(True)
plt.xlim(-2, 1)
plt.tight_layout()
plt.show()

# --------- Cálculo de parámetros ---------
ni_cm3 = 1e10  # valor típico a 300K
ni = ni_cm3 * 1e6

phi_F = Vth * np.log(N_A / ni)

V_T = V_FB + 2*phi_F + np.sqrt(4*q*eps_si*N_A*phi_F)/Cox

print("Resultados:")
print(f"Cox = {Cox:.3e} F/m^2")
print(f"phi_F = {phi_F:.3f} V")
print(f"(Asumido) V_FB = {V_FB:.3f} V")
print(f"Voltaje umbral estimado V_T = {V_T:.3f} V")
