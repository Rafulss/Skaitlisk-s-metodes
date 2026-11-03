import numpy as np
import matplotlib.pyplot as plt

N = 100
xi_max = 5

i = np.linspace(-N, N, 2 * N + 1) # Visi vesalie skaitļi no -100 lidz 100
xi = i * xi_max / N

h = xi_max / N


xi_int = xi[1:-1]
n = len(xi_int)

H = np.zeros((n, n))

f_pp = ((xi - h) - 2 * (xi) * (xi + h)) / h ** 2

print(np.size(f_pp))