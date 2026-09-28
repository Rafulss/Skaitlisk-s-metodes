import matplotlib.pyplot as plt
import numpy as np

# Vada garmus
L = 10

# Attēlums starp stabiem dots
s = L / 1.3

def f(alpha):
    return L / (2 * alpha) - np.sinh(s / (2 * alpha))

# Funkcijas atvasinājums
def df(alpha):
    return -L / (2 * alpha**2) + s * np.cosh(s / (2 * alpha)) / (2 * alpha**2)

a = 1

# Pirmais nākamais tuvinājums / lai saktu ciklu
b = a - f(a) / df(a)

# Precizitāte
delta = 1e-10

while abs(a - b) > delta:
    b = a
    a = b - f(b) / df(b)


alpha = b

print("alpha =", alpha)
print("alpha / L =", alpha / L)

