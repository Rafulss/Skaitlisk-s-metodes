import matplotlib.pyplot as plt
import numpy as np

a0 = np.array([
    np.pi/20, 2*np.pi/20, 3*np.pi/20, 4*np.pi/20, 5*np.pi/20,
    6*np.pi/20, 7*np.pi/20, 8*np.pi/20, 9*np.pi/20, 10*np.pi/20
    ])

N = 10000

q = -a0
w = a0

h = (w - q) / (N - 1)

# Izveido tukšu masīvu ar tik daudz vērtībām cik ir a0 [0. 0. 0.]
T = np.zeros(len(a0))

# a = alpha| 1/w(x)
def f(a, a0):
    return 1 / np.sqrt(np.cos(a) - np.cos(a0))

for i in range(len(a0)):
    a = np.linspace(q[i], w[i], N)

    # Neizmantojam galapunktus, jo tajos omega = 0 lai nav 1/0
    y = f(a[1:-1], a0[i])

    # Izmantoju 4.1.15:
    # h*(3/2*f2 + f3 + ... + f_(N-2) + 3/2*f_(N-1))

    # y[0] - pirmais elements, y[-1] - pēdējais elements
    T[i] = h[i] * (3/2 * y[0] + np.sum(y[1:-1]) + 3/2 * y[-1])

# lai iegūtu prasīto bezdimensionālo periodu, nav laika mērvieninības
T_t = T / (np.pi * np.sqrt(2))

print(" alpha0       T_tilde")
for i in range(len(a0)):
    #.6f cik cipari aiz kpomata
    print(f"{a0[i]:.6f}   {T_t[i]:.6f}")

plt.figure(figsize=(7, 4))
plt.plot(a0, T_t, 'o-')
plt.axhline(1, linestyle='--')
plt.xlabel(r'$\alpha_0$')
plt.ylabel(r'$\tilde{T}$')
plt.title(r'Matemātiskā svārsts')
plt.grid(True)
#plt.legend() #pie vairākām līknēm
plt.tight_layout()
#plt.savefig('', dpi=300, bbox_inches='tight')
plt.show()