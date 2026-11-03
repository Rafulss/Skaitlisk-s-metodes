import numpy as np
import matplotlib.pyplot as plt

# Dotie dati slāpekļa molekulai (N2)
# De dots viļņu skaitļos (cm^-1), pārvēršam uz SI vai izmantojam atbilstošas konstantes
De_cm = 79890.0  # cm^-1
k_angstrom = 2.6889  # A^-1
k_si = k_angstrom * 1e10  # m^-1

# Konstantes
c = 2.99792458e10  # gaismas ātrums cm/s
h = 6.62607015e-34  # J s
hbar = h / (2 * np.pi)
m_N = 14.0031 * 1.66054e-27  # viena N atoma masa kg
mu = m_N / 2.0  # reducētā masa kg

# Pārvēršam De uz Džouliem (J), lai sakārtotu mērogus
# E_h = h * c * De (ja izmanto cm^-1)
# Bet bezdimensionālajā formā mēs varam izteikt potenciālu attiecībā pret \hbar \omega
# Aprēķinām omega harmoniskajam tuvinājumam Morzes potenciālam:
# omega = k * sqrt(2 * De_dzouli / mu)
# De izteiksme enerģijas vienībās: De (J) = De_cm * h * c (precizējot vienības)
De_J = De_cm * 100 * h * c  # cm^-1 uz m^-1 un tad uz J

omega = k_si * np.sqrt(2.0 * De_J / mu)

# Bezdimensionālais parametrs Morzes potenciālam režģī:
# V_M(xi) / (hbar * omega) izteiksme:
# Tā kā x = xi * sqrt(hbar / (mu * omega)), tad kx = xi * k * sqrt(hbar / (mu * omega))
alpha_param = k_si * np.sqrt(hbar / (mu * omega))
lambda_m = np.sqrt(2 * De_J / (hbar * omega))  # Bieži sastopams Morzes parametrs

# Režģa parametri
N = 150
xi_max = 8.0
i = np.linspace(-N, N, 2 * N + 1)
xi = i * xi_max / N
d_xi = xi_max / N

xi_int = xi[1:-1]
n_int = len(xi_int)

# Morzes potenciāls bezdimensionālā formā: V_M(xi) / (hbar * omega)
# Izmantojot saistību: V_M / (hbar * omega) = (lambda_m^2 / 2) * (1 - exp(-alpha_param * xi))^2 
# vai vienkāršāk caur D_e un \hbar \omega:
V_morse_dimensionless = (De_J / (hbar * omega)) * (1.0 - np.exp(-alpha_param * xi_int))**2

# Hamiltona matrica Morzes potenciālam
H_morse = np.zeros((n_int, n_int))
diag_val = 1.0 / (d_xi ** 2) + V_morse_dimensionless
off_diag_val = -1.0 / (2.0 * d_xi ** 2)

np.fill_diagonal(H_morse, diag_val)
np.fill_diagonal(H_morse[1:], off_diag_val)
np.fill_diagonal(H_morse[:, 1:], off_diag_val)

# Diagonalizācija
eigenvalues_morse, _ = np.linalg.eigh(H_morse)

print("Pirmās 5 īpašenerģijas Morzes potenciālam ($\epsilon_n$):")
for idx in range(5):
    print(f"n = {idx}: {eigenvalues_morse[idx]:.6f}  (Harmoniskajam: {idx + 0.5})")