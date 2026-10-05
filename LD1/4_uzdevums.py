import matplotlib.pyplot as plt
import numpy as np

a0 = np.pi/3

N = np.arange(100, 100000, 100)

q = -a0
w = a0

h = (w - q) / (N - 1)

T = np.zeros(len(N))

# a = alpha| 1/w(x)
def f(a, a0):
    return 1 / np.sqrt(np.cos(a) - np.cos(a0))

# Precīzāka references vērtība no Gausa–Ležandra kvadratūras
N_ref = 100

x_ref, A_ref = np.polynomial.legendre.leggauss(N_ref)

a_ref = 2 * np.arcsin(
    np.sin(a0 / 2) * np.sin(np.pi * x_ref / 2)
)

da_dx_ref = (
    np.pi * np.sin(a0 / 2) * np.cos(np.pi * x_ref / 2)
    / np.sqrt(1 - np.sin(a0 / 2)**2 * np.sin(np.pi * x_ref / 2)**2)
)

T_ref = np.sum(A_ref * f(a_ref, a0) * da_dx_ref)
T_ref = T_ref / (np.pi * np.sqrt(2))


for i in range(len(N)):
    a = np.linspace(q, w, N[i])

    y = f(a[1:-1], a0)

    T[i] = h[i] * (3/2 * y[0] + np.sum(y[1:-1]) + 3/2 * y[-1])

T_t = T / (np.pi * np.sqrt(2))


#T_r = T_t[-1]
E = np.abs(T_t - T_ref)

lnN = np.log10(N[:-1])
lnE = np.log10(E[:-1])

#ar mazako kvadratu metodi pielago polinomu
g = np.polyfit(lnN[0:10], lnE[0:10], 1)
gamma = g[0] # taisnes slīpuma koef.
C = g[1] # brīvais b

lnE_fit = g[0] * lnN + g[1]

#print(" N       T_tilde")
#for i in range(len(N)):
#    print(f"{N[i]:.6f}   {T[i]:.6f}")

print(r'γ =', gamma, r'C =', C)

plt.figure(figsize=(7, 4))
plt.plot(lnN[:-1], lnE[:-1], 'o-', markersize=3)
plt.plot(lnN, lnE_fit, '-', linewidth=2)
#plt.axhline(1, linestyle='--')
plt.xlabel(r'$\log_{10}(N)$')
plt.ylabel(r'$\log_{10}(E_N) = |\tilde{T}_N - \tilde{T}_{ref}|$')
plt.title(r'Integrēšanas metodes konverģence')
plt.grid(True)
#plt.legend() #pie vairākām līknēm
plt.tight_layout()
#plt.savefig('', dpi=300, bbox_inches='tight')
plt.show()