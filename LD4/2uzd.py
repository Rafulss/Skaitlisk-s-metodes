import numpy as np
import matplotlib.pyplot as plt

# Parametri
L = 8
N = 200

# Rezgja solis
h = 2 * L / (N + 1)

# Koordinatas
x = np.zeros(N)

for i in range(N):
    x[i] = -L + (i + 1) * h

# Hamiltona matrica
H = np.zeros((N, N))

for i in range(N):
    H[i, i] = 2 / h**2 + x[i]**2

for i in range(N - 1):
    H[i, i + 1] = -1 / h**2
    H[i + 1, i] = -1 / h**2

# Ipasvertibas un ipasvektori
E, psi = np.linalg.eigh(H)

# Izvadam piecas mazakas ipasvertibas
print("Piecas mazakas ipasvertibas: ")

for i in range(5):
    print(E[i])

# Attelot piecas mazakas ipasvertibas
plt.plot(range(5), E[:5], "o")
plt.xlabel("Stavokla numurs")
plt.ylabel("Ipasvertiba E")
plt.grid()
plt.show()