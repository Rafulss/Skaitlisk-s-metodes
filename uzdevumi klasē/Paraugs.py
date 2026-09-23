import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return np.sqrt(x)

a = 0
b = 1
n = 3

x = np.linspace(a, b)
y = f(x)

xp = np.linspace(a, b, n)
yp = f(xp)


# Grafiks
plt.figure(figsize=(10, 5))

plt.plot(x, y)

plt.xlabel('x')
plt.ylabel('y')
plt.title(rf'')
plt.grid(True)
plt.legend()
plt.tight_layout()

#plt.savefig('', dpi=300, bbox_inches='tight')

plt.show()