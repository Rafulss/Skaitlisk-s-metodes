import matplotlib.pyplot as plt
import numpy as np

# Konstante a = (F0 * L) / (M * g) nosaka izliekumu
a = 1

N = 100

q = -2 * a
w = 2 * a

x = np.linspace(q, w, N + 1)

h = (w - q) / N

y = a * np.cosh(x / a) - a


# 3 punktu formulas 1 un 2 atvasinājums
y_1 = (y[2:] - y[:-2]) / (2 * h)
y_2 = (y[2:] - 2 * y[1:-1] + y[:-2]) / h**2

x3 = x[1:-1]

# 5 punktu formulas 1 un 2 atvasinājums
y1 = (y[:-4]- 8 * y[1:-3] + 8 * y[3:-1] - y[4:]) / (12 * h)
y2 = (-y[:-4] + 16 * y[1:-3] - 30 * y[2:-2] + 16 * y[3:-1] - y[4:]) / (12 * h**2)

x5 = x[2:-2]

j = np.sqrt(1 + (y1)**2) / a
n = np.sqrt(1 + (y_1)**2) / a

print("  x3         y_1         y_2         x5         y1         y2")

for i in range(5):
    print(f"{x3[i]:.6f}   {y_1[i]:.6f}   {y_2[i]:.6f}   {x5[i]:.6f}   {y1[i]:.6f}   {y2[i]:.6f}")

#grafiks
plt.figure(figsize=(7, 4))
# 3 punktu formula
plt.plot(x3, y_2, 'o', markersize=3, label=r"$y''$ (3 punktu formula)")
plt.plot(x3, n, '-', label=r"$\sqrt{1+(y')^2}/a$ (3 punktu formula)")

# 5 punktu formula
plt.plot(x5, y2, 's', markersize=3, label=r"$y''$ (5 punktu formula)")
plt.plot(x5, j, '--', label=r"$\sqrt{1+(y')^2}/a$ (5 punktu formula)")

plt.xlabel(r"${x}$")
plt.ylabel(r"")
plt.title(r'')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()