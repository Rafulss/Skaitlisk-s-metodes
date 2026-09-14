import matplotlib.pyplot as plt
import numpy as np

t = np.linspace(0, 1e-7, 100)

x = 1 - np.sqrt(1 - t ** 2)

plt.figure(figsize=(10, 5))
plt.plot(x)
plt.xlabel("t")
plt.ylabel("")
plt.grid()
plt.show()