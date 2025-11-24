import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# Parámetros
e_o = 88.54e-15  # F/cm
epsilon_si = 11.7 * e_o
e_ox = 3.9*e_o
chi_si = 4.05  #eV
chi_ox  =0.95 #eV
Eg = 1.12 #eV
tox = 16e-9  #en nm 
q = 1.602e-19
k = 1.381e-23
T = 300  # K
V_th = (k*T)/q

ni = 1e10   # cm^-3
Na = 2e16   # cm^-3

Cox = e_ox/(tox*100)  # paso a cm el tox

psi_B = -V_th * np.log(Na / ni)  # Para sustrato tipo p

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


plt.plot([Vg(_y) for _y in y], y,alpha = 0.85, linewidth=4.5, color = "orange",label="Tensión de Superficie vs VG")

plt.ylabel(r'$\psi_s$ (V)')
plt.xlabel(r'$V_G$ (V)')
plt.title(r"Curva de $\psi_s$ vs $V_G$")

plt.xlim(-2, 1)
plt.legend()

# plt.ylim(0, 1.2)
plt.grid()
plt.show()



#Capacitancia en funcion de V
# Vamos a recalcular Cd con la fórmula correcta para la región de interés:
# --- Gráfica ---

# Definición de tu eje X (Voltaje de Puerta)
VG_values = np.linspace(-2, 1, 2000)

#VT = VF B − 2ψB + γ raiz{−2ψB}

factor_body = np.sqrt(2*q*Na*epsilon_si)/Cox

VT = VFB - 2*psi_B + factor_body*np.sqrt(-2*psi_B)

# Aquí C_mos_total es igual a la constante Cox
C_cox_vector = np.full_like(VG_values, Cox)

print("Cuanto vale Cox = ", C_cox_vector*1e9)
print("Cuanto vale VT = ",VT)



C_mos_total_curva = np.zeros_like(VG_values) # Array para almacenar la curva C_MOS

# Capacitancia Mínima (C_min)


# 1. Acumulación (VG <= VFB)
mask_acumulacion = VG_values <= VFB
C_mos_total_curva[mask_acumulacion] = Cox

# 2. Inversión (VG >= VT) - Alta Frecuencia
mask_inversion = VG_values >= VT
C_mos_total_curva[mask_inversion] = Cox

# 3. Agotamiento (VFB < VG < VT)





# Modelo de vaciamiento
mask_dep = (VG_values > VFB) & (VG_values < VT)

argumento_vaciamiento = (4/(factor_body**2))*(VG_values - VFB) + 1 



C_mos_total_curva[mask_dep] = Cox / np.sqrt(argumento_vaciamiento[mask_dep])


plt.figure(figsize=(7,5))

# **CORRECCIÓN:** Multiplica C_mos_total por 1e6 para graficar en uF/cm^2
plt.plot(VG_values, C_cox_vector*1e9 , alpha=0.85,linestyle="--" ,linewidth=4.5, label=r'$C_{ox}$')
plt.plot(VG_values, C_mos_total_curva*1e9 , alpha=0.85 ,linewidth=4.5, color = "orange" , label=r'$C$')

plt.xlabel(r'$V_G$ (V)')
# **CORRECCIÓN:** Ajusta la etiqueta del eje Y a las nuevas unidades
plt.ylabel(r'$C$ ($ nF/cm^2$)') 
plt.title('C vs VG')

# Recta Vertical en VFB
plt.axvline(x=VFB, 
            color='gray', 
            linestyle='-.', 
            linewidth=2, 
            label=r'$V_{FB}$')

# Recta Vertical en VT
plt.axvline(x=VT, 
            color='purple', 
            linestyle='-.', 
            linewidth=2, 
            label=r'$V_T$')

# Para ver mejor la línea, ajustamos los límites del eje Y
# Cox en uF/cm^2 es aprox 2.15
plt.ylim(0, 300) 

plt.grid(True, which="both", linestyle="--")
plt.legend()
plt.tight_layout()
plt.show()



