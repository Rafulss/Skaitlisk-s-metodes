import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. kopa: N izmaiņas (xi_max = 5)
# ==========================================

# 1. variants: N = 10, xi_max = 5
N = 10         
xi_max = 5

i_vals = np.arange(-N, N + 1)
dxi = xi_max / N
xi = i_vals * dxi

dim = 2 * N + 1
H = np.zeros((dim, dim))
coeff = 1.0 / (dxi ** 2)

for i in range(dim):
    H[i, i] = coeff + 0.5 * (xi[i] ** 2)
    if i > 0:
        H[i, i - 1] = -0.5 * coeff
    if i < dim - 1:
        H[i, i + 1] = -0.5 * coeff

eigenvalues, eigenvectors = np.linalg.eigh(H)

for n in range(5):
    norm_factor = np.sqrt(np.sum(eigenvectors[:, n] ** 2) * dxi)
    eigenvectors[:, n] /= norm_factor


# 2. variants (_1): N_1 = 100, xi_max_1 = 5
N_1 = 100         
xi_max_1 = 5

i_vals_1 = np.arange(-N_1, N_1 + 1)
dxi_1 = xi_max_1 / N_1
xi_1 = i_vals_1 * dxi_1

dim_1 = 2 * N_1 + 1
H_1 = np.zeros((dim_1, dim_1))
coeff_1 = 1.0 / (dxi_1 ** 2)

for i in range(dim_1):
    H_1[i, i] = coeff_1 + 0.5 * (xi_1[i] ** 2)
    if i > 0:
        H_1[i, i - 1] = -0.5 * coeff_1
    if i < dim_1 - 1:
        H_1[i, i + 1] = -0.5 * coeff_1

eigenvalues_1, eigenvectors_1 = np.linalg.eigh(H_1)

for n in range(5):
    norm_factor_1 = np.sqrt(np.sum(eigenvectors_1[:, n] ** 2) * dxi_1)
    eigenvectors_1[:, n] /= norm_factor_1


# 3. variants (_2): N_2 = 1000, xi_max_2 = 5
N_2 = 1000         
xi_max_2 = 5

i_vals_2 = np.arange(-N_2, N_2 + 1)
dxi_2 = xi_max_2 / N_2
xi_2 = i_vals_2 * dxi_2

dim_2 = 2 * N_2 + 1
H_2 = np.zeros((dim_2, dim_2))
coeff_2 = 1.0 / (dxi_2 ** 2)

for i in range(dim_2):
    H_2[i, i] = coeff_2 + 0.5 * (xi_2[i] ** 2)
    if i > 0:
        H_2[i, i - 1] = -0.5 * coeff_2
    if i < dim_2 - 1:
        H_2[i, i + 1] = -0.5 * coeff_2

eigenvalues_2, eigenvectors_2 = np.linalg.eigh(H_2)

for n in range(5):
    norm_factor_2 = np.sqrt(np.sum(eigenvectors_2[:, n] ** 2) * dxi_2)
    eigenvectors_2[:, n] /= norm_factor_2


# ==========================================
# 2. kopa: xi_max izmaiņas (N = 100)
# ==========================================

# 4. variants (_3): N_3 = 100, xi_max_3 = 1
N_3 = 100         
xi_max_3 = 1

i_vals_3 = np.arange(-N_3, N_3 + 1)
dxi_3 = xi_max_3 / N_3
xi_3 = i_vals_3 * dxi_3

dim_3 = 2 * N_3 + 1
H_3 = np.zeros((dim_3, dim_3))
coeff_3 = 1.0 / (dxi_3 ** 2)

for i in range(dim_3):
    H_3[i, i] = coeff_3 + 0.5 * (xi_3[i] ** 2)
    if i > 0:
        H_3[i, i - 1] = -0.5 * coeff_3
    if i < dim_3 - 1:
        H_3[i, i + 1] = -0.5 * coeff_3

eigenvalues_3, eigenvectors_3 = np.linalg.eigh(H_3)

for n in range(5):
    norm_factor_3 = np.sqrt(np.sum(eigenvectors_3[:, n] ** 2) * dxi_3)
    eigenvectors_3[:, n] /= norm_factor_3


# 5. variants (_4): N_4 = 100, xi_max_4 = 10
N_4 = 100         
xi_max_4 = 10

i_vals_4 = np.arange(-N_4, N_4 + 1)
dxi_4 = xi_max_4 / N_4
xi_4 = i_vals_4 * dxi_4

dim_4 = 2 * N_4 + 1
H_4 = np.zeros((dim_4, dim_4))
coeff_4 = 1.0 / (dxi_4 ** 2)

for i in range(dim_4):
    H_4[i, i] = coeff_4 + 0.5 * (xi_4[i] ** 2)
    if i > 0:
        H_4[i, i - 1] = -0.5 * coeff_4
    if i < dim_4 - 1:
        H_4[i, i + 1] = -0.5 * coeff_4

eigenvalues_4, eigenvectors_4 = np.linalg.eigh(H_4)

for n in range(5):
    norm_factor_4 = np.sqrt(np.sum(eigenvectors_4[:, n] ** 2) * dxi_4)
    eigenvectors_4[:, n] /= norm_factor_4


# 6. variants (_5): N_5 = 100, xi_max_5 = 20
N_5 = 100         
xi_max_5 = 20

i_vals_5 = np.arange(-N_5, N_5 + 1)
dxi_5 = xi_max_5 / N_5
xi_5 = i_vals_5 * dxi_5

dim_5 = 2 * N_5 + 1
H_5 = np.zeros((dim_5, dim_5))
coeff_5 = 1.0 / (dxi_5 ** 2)

for i in range(dim_5):
    H_5[i, i] = coeff_5 + 0.5 * (xi_5[i] ** 2)
    if i > 0:
        H_5[i, i - 1] = -0.5 * coeff_5
    if i < dim_5 - 1:
        H_5[i, i + 1] = -0.5 * coeff_5

eigenvalues_5, eigenvectors_5 = np.linalg.eigh(H_5)

for n in range(5):
    norm_factor_5 = np.sqrt(np.sum(eigenvectors_5[:, n] ** 2) * dxi_5)
    eigenvectors_5[:, n] /= norm_factor_5


# ==========================================
# Attēlošana: Augsta kontrasta krāsas un biezas līnijas
# ==========================================
scale_factor = 0.4

# ------------------------------------------
# 1. Grafiks: N izmaiņas ($N = 10, 100, 1000$)
# ------------------------------------------
plt.figure(figsize=(6, 4))

for n in range(5):
    offset = eigenvalues[n]
    
    plt.plot(xi, scale_factor * eigenvectors[:, n] + offset, color='#0055ff', linewidth=2.2, label=f'N={N}' if n==0 else "")
    plt.plot(xi_1, scale_factor * eigenvectors_1[:, n] + offset, color='#ff6600', linewidth=2.2, linestyle='--', label=f'N_1={N_1}' if n==0 else "")
    plt.plot(xi_2, scale_factor * eigenvectors_2[:, n] + offset, color='#009933', linewidth=2.2, linestyle=':', label=f'N_2={N_2}' if n==0 else "")
    
    plt.axhline(y=offset, color='gray', linestyle='--', linewidth=0.7, alpha=0.7)

plt.title("N izmaiņu ietekme ($N = 10, 100, 1000$), $\\xi_{max} = 5$", fontsize=12, fontweight='bold')
plt.xlabel("Koordināta $\\xi$", fontsize=11)
plt.ylabel("Skalēta $\\psi_n(\\xi) + \\varepsilon_n$", fontsize=11)
plt.xlim(-xi_max, xi_max)
plt.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

# ------------------------------------------
# 2. Grafiks: $\\xi_{max}$ izmaiņas ($\\xi_{max} = 1, 10, 20$)
# ------------------------------------------
plt.figure(figsize=(6, 4))

for n in range(5):
    offset = eigenvalues_4[n]
    
    plt.plot(xi_3, scale_factor * eigenvectors_3[:, n] + offset, color='#9900cc', linewidth=2.2, label=f'$\\xi_{{max}}$=1' if n==0 else "")
    plt.plot(xi_4, scale_factor * eigenvectors_4[:, n] + offset, color='#cc0000', linewidth=2.2, linestyle='--', label=f'$\\xi_{{max}}$=10' if n==0 else "")
    plt.plot(xi_5, scale_factor * eigenvectors_5[:, n] + offset, color='#008b8b', linewidth=2.2, linestyle=':', label=f'$\\xi_{{max}}$=20' if n==0 else "")
    
    plt.axhline(y=offset, color='gray', linestyle='--', linewidth=0.7, alpha=0.7)

plt.title("$\\xi_{max}$ izmaiņu ietekme ($\\xi_{max} = 1, 10, 20$), $N = 100$", fontsize=12, fontweight='bold')
plt.xlabel("Koordināta $\\xi$", fontsize=11)
plt.ylabel("Skalēta $\\psi_n(\\xi) + \\varepsilon_n$", fontsize=11)
plt.xlim(-20, 20)
plt.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()