import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk

# Sākuma novirzes leņķi (radiānos)
a0 = np.arange(1, 11) * np.pi / 20

# Gausa–Ležandra kvadratūras punktu skaits
N = 32

# Integrālfunkcija
def f(a, a0):
    return 1 / np.sqrt(np.cos(a) - np.cos(a0))

# Gausa–Ležandra mezgli un svari intervālā [-1, 1]
x, A = np.polynomial.legendre.leggauss(N)

# Aprēķinām periodu katrai sākuma novirzei
T = np.zeros(len(a0))

for i in range(len(a0)):
    # Mainīgā aizvietošana:
    # sin(a/2) = sin(a0/2) * sin(pi*x/2)
    a = 2 * np.arcsin(
        np.sin(a0[i] / 2) * np.sin(np.pi * x / 2)
    )

    # da/dx
    da_dx = (
        np.pi * np.sin(a0[i] / 2) * np.cos(np.pi * x / 2)
        / np.sqrt(
            1 - np.sin(a0[i] / 2)**2 * np.sin(np.pi * x / 2)**2
        )
    )

    # Integrālis
    T[i] = np.sum(A * f(a, a0[i]) * da_dx)

# Bezdimensionālais skaitliskais periods
T_t = T / (np.pi * np.sqrt(2))

# Teorētiskā līkne:
# T_t = (2/pi) * K(sin^2(a0/2))
# ellipk funkcijai padod parametru m = k^2
a0_teor = np.linspace(0, np.pi / 2, 400)
T_t_teor = (2 / np.pi) * ellipk(np.sin(a0_teor / 2)**2)

# Rezultātu izdrukāšana
print(" alpha0 (rad)     T_t skaitliski")
for leņķis, periods in zip(a0, T_t):
    print(f" {leņķis:10.6f}     {periods:.10f}")

# Grafiks
plt.figure(figsize=(7, 4))
plt.plot(a0, T_t, "o", label="Skaitliskais rezultāts")
plt.plot(
    a0_teor,
    T_t_teor,
    "-",
    label=r"Teorija: $\tilde{T}=\frac{2}{\pi}K(\sin^2(\alpha_0/2))$"
)
plt.axhline(
    1,
    linestyle="--",
    color="gray",
    label="Mazu svārstību robeža"
)

plt.xlabel(r"$\alpha_0$")
plt.ylabel(r"$\tilde{T}$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()