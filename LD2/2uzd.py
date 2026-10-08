
# grafiks starpībai

import matplotlib.pyplot as plt
import numpy as np

# Konstante a = (F0 * L) / (M * g) nosaka izliekumu
a = 1
N = 100

x0 = -2 * a
xN = 2 * a

# Formulas elements
p = 0

# Izliekums
h = (xN - x0) / N

x = np.linspace(x0, xN, N + 1)

def y(x):
    return a * np.cosh(x / a) - a

# 3 un 5 punktu prim
def p3_yp(x, h):
    return 1 / h * ((p - 1 / 2) * y(x - h) - 2 * p * y(x) + (p + 1 / 2) * y(x + h))
def p5_yp(x, h):
    return 1 / h * (((2 * p **3 - 3*p**2 - p + 1) / 12) * y(x - 2*h) - ((4*p**3 - 3*p**2 - 8*p + 4) / 6) * y(x - h) + ((2*p**3 - 5*p) / 2) * y(x) + ((2*p**3 + 3*p**2 - p - 1) / 12) * y(x + 2*h) - ((4*p**3 + 3*p**2 - 8*p - 4) / 6) * y(x + h))

# 3 un 5 punktu primprim
def p3_ypp(x, h):
    return (y(x + h) - 2 * y(x) + y(x - h)) / h**2
def p5_ypp(x, h):
    return (-y(x + 2*h) + 16 * y(x + h) - 30 * y(x) + 16 * y(x - h) - y(x - 2*h)) / (12 * h**2) 

# Funkciju pierādīšanai
j3 = np.sqrt(1 + (p3_yp(x, h))**2) / a
j5 = np.sqrt(1 + (p5_yp(x, h))**2) / a

#grafiks
plt.figure(figsize=(7, 4))
# 3 punktu formula
plt.plot(x, p3_ypp(x, h), label="3 punktu: y''")
plt.plot(x, j3, "x", label="3 punktu: sqrt(1+(y')²)/a")

# 5 punktu formula
plt.plot(x, p5_ypp(x, h), label="5 punktu: y''")
plt.plot(x, j5, "o", markersize=3, label="5 punktu: sqrt(1+(y')²)/a")

plt.xlabel("x")
plt.ylabel("y / j primprim")
plt.title(r'')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


ypp3 = p3_ypp(x, h)
ypp5 = p5_ypp(x, h)

# 1. un 2. līnija:
# Vai katrā metodē vienādojuma kreisā puse y'' sakrīt ar labo pusi j?
starpiba_3 = np.abs(ypp3 - j3)
starpiba_5 = np.abs(ypp5 - j5)

# 3. un 4. līnija:
# Kā atšķiras 3 un 5 punktu formulas savā starpā?
starpiba_ypp = np.abs(ypp3 - ypp5)
starpiba_j = np.abs(j3 - j5)

plt.figure(figsize=(9, 5))

plt.semilogy(x, starpiba_3,
             label=r"3 punktu: $|y'' - \sqrt{1+(y')^2}/a|$")

plt.semilogy(x, starpiba_5,
             label=r"5 punktu: $|y'' - \sqrt{1+(y')^2}/a|$")

plt.semilogy(x, starpiba_ypp,
             label=r"Starpība starp 3 un 5 punktu $y''$")

plt.semilogy(x, starpiba_j,
             label=r"Starpība starp 3 un 5 punktu $\sqrt{1+(y')^2}/a$")

plt.xlabel("x")
plt.ylabel("Absolūtā starpība")
#plt.title("Starpības starp 3 un 5 punktu metodēm")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()