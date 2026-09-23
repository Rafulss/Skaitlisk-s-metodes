import matplotlib.pyplot as plt
import numpy as np

# Dotā funkcija
def f(x):
    return np.sqrt(x)

# Intervāls
a = 0
b = 1

# Punktu skaits interpolācijai (var izvēlēties no 2 līdz 5)
n = 30

# Blīvs punktu masīvs funkcijas un polinoma zīmēšanai
x = np.linspace(a, b)
y = f(x)

# Interpolācijas punkti
xp = np.linspace(a, b, n)
yp = f(xp)

# Lagranža interpolācijas polinoms
def lagrange(x, xp, yp):
    result = np.zeros_like(x, dtype=float)

    for j in range(len(xp)):
        Lj = np.ones_like(x, dtype=float)

        for i in range(len(xp)):
            if i != j:
                Lj *= (x - xp[i]) / (xp[j] - xp[i])

        result += yp[j] * Lj

    return result

# Aprēķinām interpolējošā polinoma vērtības
yl = lagrange(x, xp, yp)

# Grafiks
plt.figure(figsize=(10, 5))

plt.plot(x, y, linewidth=2, label=r'Dotā funkcija')
plt.plot(x, yl, '--', linewidth=2,
         label=f'Lagranža polinoms ({n} punkti)')
plt.plot(xp, yp, 'ro', markersize=7, label='Interpolācijas punkti')

plt.xlabel('x')
plt.ylabel('y')
plt.title(rf'Lagranža interpolācija funkcija')
plt.grid(True)
plt.legend()
plt.tight_layout()

#plt.savefig('lagranza_sin_interpolacija.png', dpi=300, bbox_inches='tight')

plt.show()