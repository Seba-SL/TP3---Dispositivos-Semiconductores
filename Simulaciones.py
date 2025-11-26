import numpy as np
import matplotlib.pyplot as plt

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

psi_bi =  (Eg/2) + V_th*np.log(Na/ni) 

VFB = - psi_bi 

factor_body = np.sqrt(2*q*Na*epsilon_si)/Cox

VT = VFB - 2*psi_B + factor_body*np.sqrt(-2*psi_B)




#Capacitancia en funcion de V
# Vamos a recalcular Cd con la fórmula correcta para la región de interés:
# --- Gráfica ---

# Definición de tu eje X (Voltaje de Puerta)
VG_values = np.linspace(-2, 1, 2000)

#VT = VF B − 2ψB + γ raiz{−2ψB}


# Aquí C_mos_total es igual a la constante Cox
C_cox_vector = np.full_like(VG_values, Cox)

print("Cuanto vale Cox = ", C_cox_vector*1e9)
print("Cuanto vale VT = ",VT)
print("Cuanto vale VFB = ",VFB)




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






# Para ver mejor la línea, ajustamos los límites del eje Y
# Cox en uF/cm^2 es aprox 2.15



# Cargar datos desde el archivo
# Usamos skiprows=1 porque la primera línea es un encabezado
data = np.loadtxt("/home/sebastian/FIUBA-2do Cuatrimestre/Dispositivos Semiconductores/TP3/Datos Spice/datos.txt", skiprows=1)

# Extraer columnas
VG = data[:, 1]                      # Columna "time" = VG
curve1 = data[:, 2]                  # (1/7.225n)*d(Ig(M1))/d(Vg)
curve2 = data[:, 3]                  # (1/7.225n)*d(Ig(M2))/d(Vg)



# Graficar
plt.figure(figsize=(8,5))
plt.plot(VG, curve1 + Cox, label= "C'1",alpha = 0.85,linewidth = 4.5)
plt.plot(VG, curve2+Cox, label="C'2",alpha = 0.85,linewidth = 4.5)
# **CORRECCIÓN:** Multiplica C_mos_total por 1e6 para graficar en uF/cm^2
plt.plot(VG_values, C_cox_vector , alpha=0.7,linestyle="--" ,linewidth=4.5, label=r"$C'_{ox}$",color = "orange")
# plt.plot(VG_values, C_mos_total_curva, alpha=0.8 ,linewidth=5, color = "blue" , label=r"$C'_g$")

# plt.xlabel(r'$V_G$ (V)')
# # **CORRECCIÓN:** Ajusta la etiqueta del eje Y a las nuevas unidades
# plt.ylabel(r" $C'$  ($ nF/cm^2$ )  ") 
# plt.title("C' vs VG")

# Recta Vertical en VFB
plt.axvline(x=VFB, 
            color='gray', 
            linestyle='-.', 
            linewidth=2.5, 
            label=r'$V_{FB}$')

# Recta Vertical en VT
plt.axvline(x=VT, 
            color='purple', 
            linestyle='-.', 
            linewidth=2.5, 
            label=r'$V_T$')
plt.xticks(
    list(plt.xticks()[0]) + [VFB, VT],
    list(map(lambda x: f"{x:.2f}", plt.xticks()[0])) + [r"$V_{FB}$", r"$V_T$"],
    fontsize=14
)


plt.xlabel("VG (V)")
plt.ylabel("C/( V m^2 ) ")
plt.title("C'(VG)")
plt.grid(True)

plt.legend(fontsize = 15)
plt.tight_layout()
plt.show()
