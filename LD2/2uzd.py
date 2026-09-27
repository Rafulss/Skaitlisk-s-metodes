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