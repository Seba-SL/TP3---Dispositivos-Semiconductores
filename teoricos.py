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
Na = 2e16   # Nch = Nsus en cm^-3

Cox = e_ox/(tox*100)  # paso a cm el tox

psi_B = -V_th * np.log(Na / ni)  # Para sustrato tipo p

psi_bi = (Eg/2) + V_th*np.log(Na/ni) 

VFB = - psi_bi 

factor_body = np.sqrt(2*q*Na*epsilon_si)/Cox

VT = VFB - 2*psi_B + factor_body*np.sqrt(-2*psi_B)


# Rango psi_s
psi_s_values = np.linspace(-0.2, 0.8, 200000)

# Expresión exacta
A = np.exp(-psi_s_values / V_th) + psi_s_values / V_th - 1
B = (ni**2 / Na**2) * (np.exp(psi_s_values / V_th) - psi_s_values / V_th - 1)

argument = A + B
argument[argument < 0] = np.nan

Q_s = np.sqrt(2 * epsilon_si * q * V_th * Na) * np.sqrt(argument)
Q_s_abs = np.abs(Q_s)
         
phi_F = -V_th * np.log(Na / ni)
psi_inv =- 2 * phi_F          





# ============================================================
#   GRAFICO DE |Q's| EN FUNCION DE psi_s 
# ============================================================

# -----------------------
# APROXIMACIONES SOLO EN SUS RANGOS DE VALIDEZ
# -----------------------

# ---- 1) Acumulación: psi_s < 0
mask_acumulacion = psi_s_values < 0
Q_s_acumulacion = np.abs(np.sqrt(2 * epsilon_si * q * V_th * Na)* (-psi_s_values/(2*V_th)))
Q_s_acumulacion[~mask_acumulacion] = np.nan

# ---- 2) Deplexión/vaciamiento: 0 < psi_s < 2 psi_B
mask_vaciamiento = (psi_s_values > 0) & (psi_s_values < -2*psi_B)
Q_s_vaciamiento = np.abs(np.sqrt(2 * epsilon_si * q * Na * psi_s_values))
Q_s_vaciamiento[~mask_vaciamiento] = np.nan

# ---- 3) Inversión fuerte: psi_s > 2 psi_B
mask_inversion = psi_s_values > -2*psi_B
Q_s_inversion = np.abs(np.sqrt(2 * epsilon_si * q * Na * psi_s_values+ 2 * epsilon_si * q * V_th * (ni**2/Na)  * np.exp(psi_s_values/V_th)))
Q_s_inversion[~mask_inversion] = np.nan


# ---- Gráfica ----
plt.figure(figsize=(6,4))
plt.plot(psi_s_values, Q_s_abs, alpha=0.7, linewidth=4.5, label="Exacto")
plt.plot(psi_s_values, Q_s_acumulacion,alpha=0.7, linewidth=3.5,linestyle="--", color="red",label="Aprox. Acumulación (ψ_s < 0)")
plt.plot(psi_s_values, Q_s_vaciamiento,alpha=0.7, linewidth=3.5, linestyle="--",color="green",label="Aprox. Vaciamiento (0 < ψ_s < 2ψ_B)")
plt.plot(psi_s_values, Q_s_inversion, alpha=0.7,linewidth=3.5,linestyle="--", color="orange",label="Aprox. Inversión (ψ_s > 2ψ_B)")


# ----------------------------------------------------
# MARCAR REGIONES 
# ----------------------------------------------------
plt.axvline(0, color="gray", linestyle=":", linewidth=2)
plt.text(0.02, 0, "Vaciamiento\nψs=0", color="gray")

plt.axvline(psi_inv, color="purple", linestyle="--", linewidth=2)
plt.text(psi_inv+0.02, 0, r"Inversión fuerte", color="purple")

plt.axvspan(psi_s_values.min(), 0, color="yellow", alpha=0.15, label="Acumulación")
plt.axvspan(0, psi_inv, color="green", alpha=0.12, label="Vaciamiento")
plt.axvspan(psi_inv, psi_s_values.max(), color="red", alpha=0.08, label="Inversión")


plt.yscale("log")
plt.xlabel(r'$\psi_s$ (V)')
plt.ylabel(r'$|Q´_s|$ (C/cm$^2$)')
plt.title(r"Curva de $|Q'_s|$ vs $\psi_s$")

plt.grid(True, which="both")
plt.legend(loc="upper right")

# Cartel
textstr = '\n'.join((
    r'$T = %d\,\mathrm{K}$' % T,
    r'$N_{sub} = N_{ch} = %.1e \mathrm{cm^{-3}}$' % Na,
    r'$\psi_B = %.3f\,\mathrm{V}$' % psi_B
))


plt.gca().text(
    0.6, 0.1,
    textstr,
    transform=plt.gca().transAxes,
    fontsize=15,
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.9)
)

plt.tight_layout()
plt.show()

















# ============================================================
#   GRAFICO DE ψ_s EN FUNCION DE VG
# ============================================================


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

plt.gca().text(0.7, 0.1, textstr, transform=plt.gca().transAxes, fontsize=15, bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))


plt.plot([Vg(_y) for _y in y], y,alpha = 0.75, linewidth=4.5, color = "blue",label="Tensión de Superficie vs VG")

# ----------------------------------------------------
# MARCAR REGIONES 
# ----------------------------------------------------
plt.axvline(0, color="gray", linestyle=":", linewidth=2)
plt.text(VFB, VT, "Vaciamiento", color="gray")
plt.axvline(VT, color="purple", linestyle="--", linewidth=2)
plt.text(VT, 1, r"Inversión fuerte", color="purple")

plt.axvspan(-2, VFB, color="yellow", alpha=0.15, label="Acumulación")
plt.axvspan(VFB, VT, color="green", alpha=0.12, label="Vaciamiento")

plt.axvspan(VT, psi_s_values.max(), color="red", alpha=0.08, label="Inversión")

plt.ylabel(r'$\psi_s$ (V)')
plt.xlabel(r'$V_G$ (V)')
plt.title(r"Curva de $\psi_s$ vs $V_G$")


plt.xticks(list(plt.xticks()[0]) + [VFB, VT],list(map(lambda x: f"{x:.2f}", plt.xticks()[0])) + [r"$V_{FB}$", r"$V_T$"],fontsize=14)



# Recta Vertical en VFB
plt.axvline(x=VFB,  color='gray', linestyle='-.', linewidth=2.5, label=f'$V_FB$={VFB:.2f} V')

# Recta Vertical en VT
plt.axvline(x=VT, 
            color='purple', 
            linestyle='-.', 
            linewidth=2.5, 
            label=f'$V_T$ = {VT:.2f} V')


plt.xlim(-2, 1)
plt.legend(fontsize=15) 
# plt.ylim(0, 1.2)
plt.grid()
plt.show()


















# ============================================================
#   GRAFICO DE C EN FUNCION DE VG
# ============================================================

# Vamos a recalcular Cd con la fórmula correcta para la región de interés:
# --- Gráfica ---
# Definición de tu eje X (Voltaje de Puerta)
VG_values = np.linspace(-2, 1, 2000)
#VT = VF B − 2ψB + γ raiz{−2ψB}

# Aquí C_mos_total es igual a la constante Cox
C_cox_vector = np.full_like(VG_values, Cox)
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
plt.plot(VG_values, C_cox_vector*1e9 , alpha=0.7,linestyle="--" ,linewidth=4.5, label=f"$C'ox = {Cox*1e9:.1f} nF/cm^2 $",color = "orange")
plt.plot(VG_values, C_mos_total_curva*1e9 , alpha=0.8 ,linewidth=5, color = "blue" , label=r"$C'_g $")


# ----------------------------------------------------
# MARCAR REGIONES 
# ----------------------------------------------------
plt.axvline(0, color="gray", linestyle=":", linewidth=2)
plt.text(VFB, VT, "Vaciamiento", color="gray")
plt.axvline(VT, color="purple", linestyle="--", linewidth=2)
plt.text(VT, 1, r"Inversión fuerte", color="purple")

plt.axvspan(-2, VFB, color="yellow", alpha=0.15, label="Acumulación")
plt.axvspan(VFB, VT, color="green", alpha=0.12, label="Vaciamiento")

plt.axvspan(VT, psi_s_values.max(), color="red", alpha=0.08, label="Inversión")

plt.xlabel(r'$V_G$ (V)')
# **CORRECCIÓN:** Ajusta la etiqueta del eje Y a las nuevas unidades
plt.ylabel(r" $C'$  ($ nF/cm^2$ )  ") 
plt.title("C' vs VG")

# Recta Vertical en VFB
plt.axvline(x=VFB, color='gray', linestyle='-.', linewidth=2.5, label=f'$VFB = {VFB:.3f} V $')

# Recta Vertical en VT
plt.axvline(x=VT, color='purple',  linestyle='-.', linewidth=2.5, label=f'$VT = {VT:.3f} V $')

plt.xticks(list(plt.xticks()[0]) + [VFB, VT],list(map(lambda x: f"{x:.1f}", plt.xticks()[0])) + [r"$V_{FB}$", r"$V_T$"], fontsize=14)

# Para ver mejor la línea, ajustamos los límites del eje Y
# Cox en uF/cm^2 es aprox 2.15


plt.ylim(0, 250) 

plt.grid(True, which="both", linestyle="--")
plt.legend(fontsize=15)  
plt.tight_layout()
plt.show()















# ============================================================
#   GRAFICO DE Qd Y Qi EN FUNCION DE psi_s  (psi_s > 0)
# ============================================================

Qi = lambda psi_s: - Cox * (Vg(psi_s) - (VFB + psi_s + factor_body * np.sqrt(psi_s)))
Qd = lambda psi_s: - np.sqrt(2 * epsilon_si * q * Na * psi_s)
# Rango psi_s
psi_s_values_QiQd = np.linspace(0, 0.8, 200000)
Qi_vals = np.array([-Qi(x) for x in psi_s_values_QiQd])
Qd_vals = np.array([-Qd(x) for x in psi_s_values_QiQd])
cruce_Q = 0.479

# ------------------------------------------------------------
#   GRAFICO
# ------------------------------------------------------------

plt.plot(psi_s_values_QiQd, Qi_vals, c='b', label=r"$|Q'_i|$", alpha=0.8, linewidth=4.5)
plt.plot(psi_s_values_QiQd, Qd_vals, c='r', label=r"$|Q'_d|$", alpha=0.8, linewidth=4.5)
plt.plot(psi_s_values_QiQd, Qi_vals + Qd_vals, c='g', label=r"$|Q'_s|$", alpha=0.8, linewidth=4.5)

plt.axvspan(cruce_Q, psi_s_values_QiQd.max(),color='violet', alpha=0.35, label=r"Región: $|Q´_i |$ > $|Q´_d |$")

# Recta Vertical en VFB
plt.axvline(x=cruce_Q, color='gray',  linestyle='-.', linewidth=2.5,  label=f"Cruce : $ |Q´_i| \;  = \; |Q´_d|$ , $\Psi_s $ = {cruce_Q} V")

plt.yscale('log')
plt.xlabel(r"$\psi_s$  [V]")
plt.ylabel(r"Densidad de carga  [C/cm$^2$]")
plt.legend(fontsize=15)
plt.grid()
plt.show()
