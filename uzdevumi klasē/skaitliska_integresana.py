import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return np.sin(x)

a = 0
b = np.pi / 2

x = np.linspace(a, b)

# Trapeču formula
I_trapece = (b - a) / 2 * (f(a) + f(b))

# Simpsona formula
x_vidus = (a + b) / 2

y0 = f(a)
y1 = f(x_vidus)
y2 = f(b)

delta = x_vidus - a

I_simpsons = (y0 / 3 + 4 * y1 / 3 + y2 / 3) * delta

# Simpsona 3/8 formula
h_38 = (b - a) / 3
x1 = a + h_38
x2 = a + 2 * h_38

I_simpsons_38 = 3 * h_38 / 8 * (
    f(a) + 3 * f(x1) + 3 * f(x2) + f(b)
)

#N apakšintervāliem
N = 10
x = np.linspace(a, b, N + 1)
y = f(x)
h = (b - a) / N

I_trapece_N = h * (
    0.5 * y[0] + np.sum(y[1:-1]) + 0.5 * y[-1]
)

print(f"Trapeču formula:                 {I_trapece:.10f}")
print(f"Simpsona 1/3 formula:            {I_simpsons:.10f}")
print(f"Simpsona 3/8 formula:            {I_simpsons_38:.10f}")
print(f"Trapeču formula, N = {N}:          {I_trapece_N:.10f}")

plt.figure(figsize=(10, 5))
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title(rf'')
plt.grid(True)
#plt.legend()
plt.tight_layout()

#plt.savefig('', dpi=300, bbox_inches='tight')

#plt.show()