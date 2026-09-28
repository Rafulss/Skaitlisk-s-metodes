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
a = 1
b = a - f(a) / df(a)

while abs(a - b) > delta:
    a = b
    b = a - f(a) / df(a)

alpha_precizs = b


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
a = 1
nut_kluda = []

b = a - f(a) / df(a)

while abs(a - b) > delta:
    nut_kluda.append(abs(b - alpha_precizs))
    a = b
    b = a - f(a) / df(a)


# Sekanšu metode
a = 2
b = 4
sek_kluda = []

while abs(b - a) > delta:
    c = (f(b) * a - f(a) * b) / (f(b) - f(a))

    a = b
    b = c

    sek_kluda.append(abs(b - alpha_precizs))


# Grafiks
plt.plot(bis_kluda, label="Bisekcija")
plt.plot(nut_kluda, label="Ņūtons")
plt.plot(sek_kluda, label="Sekantes")

plt.yscale("log")
plt.xlabel("Iterācija")
plt.ylabel("Kļūda")
plt.grid()
plt.legend()
plt.show()

print("Precīzais alpha =", alpha_precizs)