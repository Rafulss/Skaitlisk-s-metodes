import matplotlib.pyplot as plt
import numpy as np

a0 = np.pi / 3

# N tagad ir Gausa–Ležandra punktu skaits
N = np.arange(2, 21)

q = -a0
w = a0

T = np.zeros(len(N))

# a = alpha | 1 / omega(alpha)
def f(a, a0):
    return 1 / np.sqrt(np.cos(a) - np.cos(a0))

for i in range(len(N)):

    # Gausa–Ležandra mezgli un svari
    x, A = np.polynomial.legendre.leggauss(N[i])

    # Mainīgā aizvietošana:
    # sin(a/2) = sin(a0/2) * sin(pi*x/2)
    a = 2 * np.arcsin(np.sin(a0 / 2) * np.sin(np.pi * x / 2))

    # da/dx
    da_dx = (
        np.pi * np.sin(a0 / 2) * np.cos(np.pi * x / 2)
        / np.sqrt(1 - np.sin(a0 / 2)**2 * np.sin(np.pi * x / 2)**2)
    )

    # Gausa–Ležandra kvadratūra
    T[i] = np.sum(A * f(a, a0) * da_dx)

T_t = T / (np.pi * np.sqrt(2))

# References vērtība ar lielāku Gausa punktu skaitu
N_ref = 100

x_ref, A_ref = np.polynomial.legendre.leggauss(N_ref)

a_ref = 2 * np.arcsin(np.sin(a0 / 2) * np.sin(np.pi * x_ref / 2))

da_dx_ref = (
    np.pi * np.sin(a0 / 2) * np.cos(np.pi * x_ref / 2)
    / np.sqrt(1 - np.sin(a0 / 2)**2 * np.sin(np.pi * x_ref / 2)**2)
)

T_ref = np.sum(A_ref * f(a_ref, a0) * da_dx_ref)
T_ref = T_ref / (np.pi * np.sqrt(2))

E = np.abs(T_t - T_ref)

N_fit = N[:-4]
E_fit = E[:-4]

#lnN = np.log10(N_fit)
lnE = np.log10(E_fit)

# Taisnes pielāgošana pirmajiem punktiem
#g = np.polyfit(lnN, lnE, 1)
g = np.polyfit(N_fit, lnE, 1)

gamma = g[0]
C = g[1]

lnE_fit = gamma * N_fit + C

print("N       T_tilde            E")
for i in range(len(N)):
    print(f"{N[i]:2d}   {T_t[i]:.12f}   {E[i]:.3e}")

print(r'γ =', gamma)
print(r'C =', C)

plt.figure(figsize=(7, 4))
plt.plot(N_fit, lnE, 'o-', markersize=4)
plt.plot(N_fit, lnE_fit, '-', linewidth=2)

plt.xlabel(r'$N$')
plt.ylabel(r'$\log_{10}(E_N)$')
#plt.title(r'Gausa–Ležandra kvadratūras konverģence')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()