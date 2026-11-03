import numpy as np
import matplotlib.pyplot as plt

# 1. Parametru definēšana
N = 100          # punktu skaits no -N līdz N
xi_max = 10    # Maksimālā robeža koordinātai xi

# 2. Režģa izveide
i_vals = np.arange(-N, N + 1)
dxi = xi_max / N
xi = i_vals * dxi

# 3. Matricas H izveide un aizpildīšana (manuāli ar ciklu)
dim = 2 * N + 1
H = np.zeros((dim, dim))
coeff = 1.0 / (dxi ** 2)

for i in range(dim):
    # Galvenās diagonāles elements
    H[i, i] = coeff + 0.5 * (xi[i] ** 2)
    
    # Augšējā blakusdiagonāle
    if i > 0:
        H[i, i - 1] = -0.5 * coeff
        
    # Apakšējā blakusdiagonāle
    if i < dim - 1:
        H[i, i + 1] = -0.5 * coeff

# 4. Matricas diagonalizācija (īpašvērtības un īpašvektori)
eigenvalues, eigenvectors = np.linalg.eigh(H)

# Izvadām 5 mazākās īpašvērtības (enerģijas)
print("5 vismazākās īpašvērtības:", eigenvalues[:5])