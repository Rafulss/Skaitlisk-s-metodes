import numpy as np
import matplotlib.pyplot as plt

# Funkcija, kas atgriež pirmās 5 īpašvērtības dotajam N un xi_max
def calculate_eigenvalues(N, xi_max):
    i = np.linspace(-N, N, 2 * N + 1)
    xi = i * xi_max / N
    d_xi = xi_max / N
    
    xi_int = xi[1:-1]
    n_int = len(xi_int)
    
    H = np.zeros((n_int, n_int))
    diag_val = 1.0 / (d_xi ** 2) + 0.5 * (xi_int ** 2)
    off_diag_val = -1.0 / (2.0 * d_xi ** 2)
    
    np.fill_diagonal(H, diag_val)
    np.fill_diagonal(H[1:], off_diag_val)
    np.fill_diagonal(H[:, 1:], off_diag_val)
    
    eigenvalues, _ = np.linalg.eigh(H)
    return eigenvalues[:5]

# --- 1. Pētām atkarību no parametra xi_max (fiksēts N = 100) ---
N_fixed = 100
xi_max_values = np.linspace(2.0, 10.0, 30)
energies_vs_ximax = np.zeros((len(xi_max_values), 5))

for idx, xmax in enumerate(xi_max_values):
    energies_vs_ximax[idx, :] = calculate_eigenvalues(N_fixed, xmax)

# Teorētiskās vērtības harmoniskajam oscilatoram: epsilon_n = n + 0.5 (n = 0, 1, 2, 3, 4)
# Teorētiskās: 0.5, 1.5, 2.5, 3.5, 4.5

plt.figure(figsize=(10, 5))
for n in range(5):
    plt.plot(xi_max_values, energies_vs_ximax[:, n], label=f'$\\epsilon_{n}$ (skaitliski)')
plt.axhline(0.5, color='gray', linestyle='--', alpha=0.5)
plt.axhline(4.5, color='gray', linestyle='--', alpha=0.5, label='Teorētiskās robežas')
plt.title(f'Enerģiju atkarība no $\\xi_{max}$ (pie fiksēta $N = {N_fixed}$)')
plt.xlabel('$\\xi_{max}$')
plt.ylabel('Īpašenerģija $\\epsilon$')
plt.legend()
plt.grid(True)
plt.show()

# --- 2. Pētām atkarību no parametra N (fiksēts xi_max = 6.0) ---
xi_max_fixed = 6.0
N_values = np.arange(20, 200, 10)
energies_vs_N = np.zeros((len(N_values), 5))

for idx, n_pts in enumerate(N_values):
    energies_vs_N[idx, :] = calculate_eigenvalues(n_pts, xi_max_fixed)

plt.figure(figsize=(10, 5))
for n in range(5):
    plt.plot(N_values, energies_vs_N[:, n], marker='o', markersize=3, label=f'$\\epsilon_{n}$')
plt.title(f'Enerģiju konverģence atkarībā no $N$ (pie fiksēta $\\xi_{max} = {xi_max_fixed}$)')
plt.xlabel('Režģa parametrs $N$')
plt.ylabel('Īpašenerģija $\\epsilon$')
plt.legend()
plt.grid(True)
plt.show()