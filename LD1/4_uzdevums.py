import matplotlib.pyplot as plt
import numpy as np

a0 = np.pi/3

N = np.arange(100, 10100, 100)

q = -a0
w = a0

h = (w - q) / (N - 1)

T = np.zeros(len(N))

# a = alpha| 1/w(x)
def f(a, a0):
    return 1 / np.sqrt(np.cos(a) - np.cos(a0))

for i in range(len(N)):
    a = np.linspace(q, w, N[i])

    y = f(a[1:-1], a0)

    T[i] = h[i] * (3/2 * y[0] + np.sum(y[1:-1]) + 3/2 * y[-1])

T_t = T / np.pi * np.sqrt(2)


T_r = T_t[-1]
E = np.abs(T_t - T_r)

lnN = np.log(N[:-1])
lnE = np.log(E[:-1])

#ar mazako kvadratu metodi pielago polinomu
g = np.polyfit(lnN, lnE, 1)
gamma = -g[0] # taisnes slīpuma koef.
C = g[1] # brīvais b


print(" N       T_tilde")
for i in range(len(N)):
    print(f"{N[i]:.6f}   {T[i]:.6f}")

print(r'γ =', gamma, r'C =', C)

plt.figure(figsize=(7, 4))
plt.plot(N, T_t, 'o-', markersize=3)
#plt.axhline(1, linestyle='--')
plt.xlabel(r'$N$')
plt.ylabel(r'$\tilde{T}$')
plt.title(r'Bezdimensionāla metodes konverģence')
plt.grid(True)
#plt.legend() #pie vairākām līknēm
plt.tight_layout()
#plt.savefig('', dpi=300, bbox_inches='tight')
plt.show()