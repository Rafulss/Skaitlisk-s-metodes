import numpy as np
import matplotlib.pyplot as plt

# 1. Parametru definēšana
N = 100       
xi_max = 10

# 2. Režģa izveide
i_vals = np.arange(-N, N + 1)
dxi = xi_max / N
xi = i_vals * dxi

# 3. Matricas H izveide un aizpildīšana
dim = 2 * N + 1
H = np.zeros((dim, dim))
coeff = 1.0 / (dxi ** 2)

for i in range(dim):
    H[i, i] = coeff + 0.5 * (xi[i] ** 2)
    if i > 0:
        H[i, i - 1] = -0.5 * coeff
    if i < dim - 1:
        H[i, i + 1] = -0.5 * coeff

# 4. Matricas diagonalizācija
eigenvalues, eigenvectors = np.linalg.eigh(H)

# 5. Viļņu funkciju normēšana
for n in range(5):
    norm_factor = np.sqrt(np.sum(eigenvectors[:, n] ** 2) * dxi)
    eigenvectors[:, n] /= norm_factor

# 6. 4. uzdevums: Viļņu funkciju zīmēšana ar skalēšanu, lai tās nepārklātos
plt.figure(figsize=(8, 9))

scale_factor = 0.4  # Skalēšanas koeficients, kas samazina viļņu funkcijas vertikālo izmēru

for n in range(5):
    offset = eigenvalues[n]
    # Reizinām viļņu funkciju ar scale_factor un pieskaitām enerģijas līmeni (offset)
    plt.plot(xi, scale_factor * eigenvectors[:, n] + offset, label=f'n = {n}, $\\varepsilon$ = {eigenvalues[n]:.3f}')
    
    # Bāzes līnija attiecīgajai enerģijai
    plt.axhline(y=offset, color='gray', linestyle='--', linewidth=0.7)

plt.title("Harmoniskā oscilatora viļņu funkcijas")
plt.xlabel("Koordināta $\\xi$")
plt.ylabel("Skalēta $\\psi_n(\\xi) + \\varepsilon_n$")
plt.xlim(-xi_max, xi_max)
plt.legend(loc='upper right')
plt.grid(True)
plt.show()