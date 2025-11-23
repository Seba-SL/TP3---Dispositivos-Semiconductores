import numpy as np
import matplotlib.pyplot as plt

# Parámetros
e_o = 88.54e-15  # F/cm
epsilon_si = 11.7 * e_o
e_ox = 3.9*e_o



chi_si = 4.05  #eV
chi_ox  =0.95 #eV

Eg = 1.12 #eV


tox = 16e-9

q = 1.602e-19
k = 1.381e-23
T = 300  # K

V_th = (k*T)/q

ni = 1e10   # cm^-3
Na = 2e16   # cm^-3

Cox = e_ox/tox

psi_B = V_th * np.log(Na / ni)  # Para sustrato tipo p

psi_bi = chi_si*q + (Eg/2*q) + V_th*np.log(Na/ni) - chi_ox*q

VFB = - psi_bi 



# Rango psi_s
psi_s_values = np.linspace(-0.2, 0.8, 200000)

# Expresión exacta
A = np.exp(-psi_s_values / V_th) + psi_s_values / V_th - 1
B = (ni**2 / Na**2) * (np.exp(psi_s_values / V_th) - psi_s_values / V_th - 1)

argument = A + B
argument[argument < 0] = np.nan

Q_s = np.sqrt(2 * epsilon_si * q * V_th * Na) * np.sqrt(argument)
Q_s_abs = np.abs(Q_s)

# -----------------------
# APROXIMACIONES SOLO EN SUS RANGOS DE VALIDEZ
# -----------------------

# ---- 1) Acumulación: psi_s < 0
mask_acumulacion = psi_s_values < 0
Q_s_acumulacion = np.abs(np.sqrt(2 * epsilon_si * q * V_th * Na)
                         * (-psi_s_values/(2*V_th)))
Q_s_acumulacion[~mask_acumulacion] = np.nan

# ---- 2) Deplexión/vaciamiento: 0 < psi_s < 2 psi_B
mask_vaciamiento = (psi_s_values > 0) & (psi_s_values < 2*psi_B)
Q_s_vaciamiento = np.abs(np.sqrt(2 * epsilon_si * q * Na * psi_s_values))
Q_s_vaciamiento[~mask_vaciamiento] = np.nan

# ---- 3) Inversión fuerte: psi_s > 2 psi_B
mask_inversion = psi_s_values > 2*psi_B
Q_s_inversion = np.abs(np.sqrt(2 * epsilon_si * q * Na * psi_s_values
                       + 2 * epsilon_si * q * V_th * (ni**2/Na)
                       * np.exp(psi_s_values/V_th)))
Q_s_inversion[~mask_inversion] = np.nan



# ---- Gráfica ----
plt.figure(figsize=(6,4))
plt.plot(psi_s_values, Q_s_abs, alpha=0.7, linewidth=4.5, label="Exacto")
plt.plot(psi_s_values, Q_s_acumulacion,alpha=0.7, linewidth=3.5,linestyle="--", color="red",label="Aprox. Acumulación (ψ_s < 0)")
plt.plot(psi_s_values, Q_s_vaciamiento,alpha=0.7, linewidth=3.5, linestyle="--",color="green",label="Aprox. Vaciamiento (0 < ψ_s < 2ψ_B)")
plt.plot(psi_s_values, Q_s_inversion, alpha=0.7,linewidth=3.5,linestyle="--", color="orange",label="Aprox. Inversión (ψ_s > 2ψ_B)")

plt.yscale("log")
plt.xlabel(r'$\psi_s$ (V)')
plt.ylabel(r'$|Q_s|$ (C/cm$^2$)')
plt.title(r"Curva de $|Q_s|$ vs $\psi_s$")

plt.grid(True, which="both")
plt.legend()

# Cartel
textstr = '\n'.join((
    r'$T = %d\,\mathrm{K}$' % T,
    r'$N_{sub} = N_{ch} = %.1e \mathrm{cm^{-3}}$' % Na,
    r'$\psi_B = %.3f\,\mathrm{V}$' % psi_B
))

plt.gca().text(
    0.60, 0.30,
    textstr,
    transform=plt.gca().transAxes,
    fontsize=10,
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.9)
)

plt.tight_layout()
plt.show()





#psi_s

C = np.sqrt(2*epsilon_si*q*V_th*Na)
Qs = lambda phi_s: C * np.sqrt((np.e**(-phi_s / V_th) + (phi_s / V_th) - 1) + (1.01e10/Na) * (np.e**(phi_s / V_th) - (phi_s / V_th) - 1))
Vg = lambda phi_s: phi_s + (Qs(phi_s) / Cox) + VFB if phi_s > 0 else phi_s - (Qs(phi_s) / Cox) + VFB 



plt.figure(figsize=(6,4))
y = np.linspace(-2, 1, 10000)


# Cartel
textstr = '\n'.join((
    r'$T = %d\,\mathrm{K}$' % T,
    r'$N_{sub} = N_{ch} = %.1e \mathrm{cm^{-3}}$' % Na,
    r'$\psi_B = %.3f\,\mathrm{V}$' % psi_B
))

plt.gca().text(
    0.60, 0.30,
    textstr,
    transform=plt.gca().transAxes,
    fontsize=10,
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.9)
)

plt.ylabel(r'$\psi_s$ (V)')
plt.xlabel(r'$V_G$ (V)')
plt.title(r"Curva de $\psi_s$ vs $V_G$")
plt.plot([Vg(_y) for _y in y], y,alpha = 0.85, linewidth=4.5, color = "orange")
plt.xlim(-2, 1)
# plt.ylim(0, 1.2)
plt.grid()
plt.show()