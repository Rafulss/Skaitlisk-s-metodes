import matplotlib.pyplot as plt
import numpy as np

# Vada garums un attālums
L = 10
s = L / 1.3

# Precizitāte
delta = 1e-10

def f(alpha):
    return L / (2 * alpha) - np.sinh(s / (2 * alpha))

def df(alpha):
    return -L / (2 * alpha**2) + s * np.cosh(s / (2 * alpha)) / (2 * alpha**2)

# Precīzs risinājums ar Ņūtona metodi
b = 1
a = b - f(b) / df(b)

while abs(a - b) > delta:
    b = a
    a = b - f(b) / df(b)

alpha_precizs = a


# Bisekcijas metode
a = 0.1
b = 10
bis_kluda = []

while b - a > delta:
    c = (a + b) / 2
    bis_kluda.append(abs(c - alpha_precizs))

    if f(a) * f(c) < 0:
        b = c
    else:
        a = c


# Ņūtona metode
b = 1
nut_kluda = []

a = b - f(b) / df(b)

while abs(a - b) > delta:
    b = a
    a = b - f(b) / df(b)
    nut_kluda.append(abs(a - alpha_precizs))


# Sekanšu metode
a = 2
b = 3
sek_kluda = []

while abs(b - a) > delta:
    c = (f(b) * a - f(a) * b) / (f(b) - f(a))

    a = b
    b = c

    sek_kluda.append(abs(b - alpha_precizs))


# Grafiks
plt.plot(bis_kluda, "o-", label="Bisekcija")
plt.plot(nut_kluda, "o-", label="Ņūtons")
plt.plot(sek_kluda, "o-", label="Sekantes")

plt.yscale("log")
plt.xlabel("Iterācija")
plt.ylabel("Kļūda")
plt.grid()
plt.legend()
plt.show()

print("Precīzais alpha =", alpha_precizs)