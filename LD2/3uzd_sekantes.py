import matplotlib.pyplot as plt
import numpy as np

# Vada garmus
L = 10

# Attēlums starp stabiem dots
s = L / 1.3

def f(alpha):
    return L / (2 * alpha) - np.sinh(s / (2 * alpha))

# Izvēlētas a un b vērtības viena negatīva otra pozitīva
# un lai nav viena daudz lielaka par otru
a = 2
b = 4

# Precizitāte
delta = 1e-10

while abs(b - a) > delta:
    c = (f(b) * a - f(a) * b) / (f(b) - f(a))

    a = b
    b = c

alpha = b

# Pārbaudam saknes
print(f(2)) # negatīvs
print(f(4)) # pozitīvs

print("alpha =", alpha)
print("alpha / L =", alpha / L)

