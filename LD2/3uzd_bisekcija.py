import matplotlib.pyplot as plt
import numpy as np

# Vada garmus
L = 10

# Attēlums starp stabiem dots
s = L / 1.3

def f(alpha):
    return L / (2 * alpha) - np.sinh(s / (2 * alpha))

#Izvēlētas a un b vērtības
a = 0.1
b = 10

# Precizitāte
delta = 1e-10

while b - a > delta:
    c = 1 / 2 * (a + b)

    if f(a) * f(c) < 0:
        b = c
    else:
        a = c

alpha = (a + b) / 2

print("alpha =", alpha)
print("alpha / L =", alpha / L)

