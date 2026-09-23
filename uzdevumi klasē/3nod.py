import matplotlib.pyplot as plt
import numpy as np

def f(x):
    np.sin(x)

x = 0
a = -1
b = 1
c = 1/2 * (a+b)
delta = 1

if f(a) * f (c) < 0:
    b = c
else:
    a = c



# Grafiks
plt.figure(figsize=(10, 5))
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title(rf'')
plt.grid(True)
plt.tight_layout()
plt.show()