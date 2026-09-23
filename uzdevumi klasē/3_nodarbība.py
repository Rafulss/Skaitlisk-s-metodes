import matplotlib.pyplot as plt
import numpy as np

def sin_funkcija(x):
    return np.sin(x)

def cos_funkcija(x):
    return np.cos(x)

x0 = 0

precizs_sin_pirmais_atvasinajums = 1.0
precizs_cos_otrais_atvasinajums = -1.0

delta = 1 / 2**np.arange(0, 30)


sin_pirmais_3_punkti = np.zeros(len(delta))
sin_pirmais_5_punkti = np.zeros(len(delta))


cos_otrais_3_punkti = np.zeros(len(delta))
cos_otrais_5_punkti = np.zeros(len(delta))

for i in range(len(delta)):
    h = delta[i]

    sin_pirmais_3_punkti[i] = (
        sin_funkcija(x0 + h) - sin_funkcija(x0 - h)
    ) / (2 * h)

    sin_pirmais_5_punkti[i] = (
        sin_funkcija(x0 - 2*h)
        - 8 * sin_funkcija(x0 - h)
        + 8 * sin_funkcija(x0 + h)
        - sin_funkcija(x0 + 2*h)
    ) / (12 * h)


    cos_otrais_3_punkti[i] = (
        cos_funkcija(x0 + h)
        - 2 * cos_funkcija(x0)
        + cos_funkcija(x0 - h)
    ) / h**2

    cos_otrais_5_punkti[i] = (
        -cos_funkcija(x0 + 2*h)
        + 16 * cos_funkcija(x0 + h)
        - 30 * cos_funkcija(x0)
        + 16 * cos_funkcija(x0 - h)
        - cos_funkcija(x0 - 2*h)
    ) / (12 * h**2)

sin_kluda_3_punkti = np.abs(
    sin_pirmais_3_punkti - precizs_sin_pirmais_atvasinajums
)

sin_kluda_5_punkti = np.abs(
    sin_pirmais_5_punkti - precizs_sin_pirmais_atvasinajums
)

cos_kluda_3_punkti = np.abs(
    cos_otrais_3_punkti - precizs_cos_otrais_atvasinajums
)

cos_kluda_5_punkti = np.abs(
    cos_otrais_5_punkti - precizs_cos_otrais_atvasinajums
)

print("f(x) = sin(x), pirmais atvasinājums punktā x = 0")
print("Precīzā vērtība: 1")
print()
print(" delta         3 punktu formula       5 punktu formula")

for i in range(len(delta)):
    print(
        f"{delta[i]:.8e}   "
        f"{sin_pirmais_3_punkti[i]:.12f}       "
        f"{sin_pirmais_5_punkti[i]:.12f}"
    )

print()
print("f(x) = cos(x), otrais atvasinājums punktā x = 0")
print("Precīzā vērtība: -1")
print()
print(" delta         3 punktu formula       5 punktu formula")

for i in range(len(delta)):
    print(
        f"{delta[i]:.8e}   "
        f"{cos_otrais_3_punkti[i]:.12f}       "
        f"{cos_otrais_5_punkti[i]:.12f}"
    )

plt.figure(figsize=(8, 5))

# sin(x) pirmā atvasinājuma kļūdas
plt.loglog(
    delta, sin_kluda_3_punkti,
    'o-', label=r'$\sin(x)$, 1. atvasinājums, 3 punkti'
)

plt.loglog(
    delta, sin_kluda_5_punkti,
    's-', label=r'$\sin(x)$, 1. atvasinājums, 5 punkti'
)

# cos(x) otrā atvasinājuma kļūdas
plt.loglog(
    delta, cos_kluda_3_punkti,
    '^-', label=r'$\cos(x)$, 2. atvasinājums, 3 punkti'
)

plt.loglog(
    delta, cos_kluda_5_punkti,
    'd-', label=r'$\cos(x)$, 2. atvasinājums, 5 punkti'
)

plt.xlabel(r'$\delta = h$')
plt.ylabel('Absolūtā kļūda')
plt.grid(True, which='both')
plt.legend()
plt.tight_layout()
plt.show()