import matplotlib.pyplot as plt
import numpy as np

a0 = np.array([
    np.pi/20, 2*np.pi/20, 3*np.pi/20, 4*np.pi/20, 5*np.pi/20,
    6*np.pi/20, 7*np.pi/20, 8*np.pi/20, 9*np.pi/20, 10*np.pi/20
])

# Gausa–Ležandra kvadratūras punktu skaits
N = 32

q = -a0
w = a0

T = np.zeros(len(a0))

# a = alpha | 1 / omega(alpha)
def f(a, a0):
    return 1 / np.sqrt(np.cos(a) - np.cos(a0))

# Gausa–Ležandra mezgli x un svari A intervālā [-1, 1]
x, A = np.polynomial.legendre.leggauss(N)

for i in range(len(a0)):

    # Mainīgā aizvietošana:
    # sin(a/2) = sin(a0/2) * sin(pi*x/2)
    a = 2 * np.arcsin(np.sin(a0[i] / 2) * np.sin(np.pi * x / 2))

    # da/dx
    da_dx = (
        np.pi * np.sin(a0[i] / 2) * np.cos(np.pi * x / 2)
        / np.sqrt(1 - np.sin(a0[i] / 2)**2 * np.sin(np.pi * x / 2)**2)
    )

    # Integrālis: summa A_i * f(a_i, a0) * da/dx
    T[i] = np.sum(A * f(a, a0[i]) * da_dx)

# Bezdimensionālais periods
T_t = T / (np.pi * np.sqrt(2))

print(" alpha0       T_tilde")
for i in range(len(a0)):
    print(f"{a0[i]:.6f}   {T_t[i]:.10f}")

plt.figure(figsize=(7, 4))
plt.plot(a0, T_t, 'o-') 
plt.axhline(1, linestyle='--', color='black')
plt.xlabel(r'$\alpha_0$')
plt.ylabel(r'$\tilde{T}$')
#plt.title(r'Matemātiskais svārsts')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()