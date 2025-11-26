import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------
#  Variables
# ----------------------------------------------------
q     = 1.6e-19               
eps0  = 8.85e-14              
eps_s = 11.7 * eps0
eps_ox = 3.9 * eps0
Vth = 25.9e-3
N_A   = 2e16                  
ni    = 1e10                  
t_ox  = 5e-7                  

Cox = eps_ox / t_ox

Vg_min = -2
Vg_max =  1

kT_q = 0.0259                
phi_F = -kT_q * np.log(N_A / ni)
psi_inv =- 2 * phi_F          

Vfb = -phi_F
gamma = np.sqrt(2 * q * eps_s * N_A) / Cox

# ----------------------------------------------------
# Vector de simulación
psi_s = np.linspace(-0.4, 1.2, 6000)

plt.figure(figsize=(11,6))

# Colores
color_Qd = "black"
color_Qi = "tab:blue"
color_Qs = "tab:red"

# ----------------------------------------------------
# GRAFICAR SOLO Vg_min y Vg_max
# ----------------------------------------------------
for Vg, label_tag in [(Vg_min, "Vg_min"), (Vg_max, "Vg_max")]:

    # -------- Qd: solo para psi_s >= 0 --------
    Qd = np.zeros_like(psi_s)
    mask_Qd = psi_s >= 0
    Qd[mask_Qd] = -np.sqrt(2 * eps_s * q * N_A * psi_s[mask_Qd])

    # -------- Qi: solo para inversion fuerte ψs ≥ ψ_inv --------
    Qi = np.zeros_like(psi_s)
    mask_inv = psi_s >= psi_inv
    #Qi[mask_inv] = -Cox * (Vg - (Vfb + psi_s[mask_inv] + gamma*np.sqrt(np.maximum(psi_s[mask_inv], 0))))
    Qi_term = Vg - (Vfb + psi_s[mask_inv] + gamma*np.sqrt(np.maximum(psi_s[mask_inv], 0)))
    Qi[mask_inv] = -Cox * np.maximum(Qi_term, 0)   # IMPORTANTE


    # -------- Total Qs --------
    Qs = Qd + Qi

    # -------- Graficar --------
    plt.plot(psi_s, Qd, '-',  linewidth=2, color=color_Qd,
             label=f"Qd ({label_tag})")

    plt.plot(psi_s, Qi, '--', linewidth=2, color=color_Qi,
             label=f"Qi ({label_tag})")

    plt.plot(psi_s, Qs, '-',  linewidth=2, color=color_Qs,
             label=f"Qs ({label_tag})")

# ----------------------------------------------------
# MARCAR REGIONES 
# ----------------------------------------------------
plt.axvline(0, color="gray", linestyle=":", linewidth=2)
plt.text(0.02, 0, "Agotamiento\nψs=0", color="gray")

plt.axvline(psi_inv, color="purple", linestyle="--", linewidth=2)
plt.text(psi_inv+0.02, 0, r"Inversión fuerte", color="purple")

plt.axvspan(psi_s.min(), 0, color="yellow", alpha=0.15, label="Acumulación")
plt.axvspan(0, psi_inv, color="green", alpha=0.12, label="Agotamiento")
plt.axvspan(psi_inv, psi_s.max(), color="red", alpha=0.08, label="Inversión")

# ----------------------------------------------------
# Gráfico final
# ----------------------------------------------------
plt.xlabel("ψs (V)")
plt.ylabel("Densidad superficial de carga (C/cm²)")
plt.title("Qs, Qd y Qi para Vg_min y Vg_max")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
