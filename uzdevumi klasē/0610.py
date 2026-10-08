import numpy as np
import matplotlib.pyplot as plt

N = 100

x = np.linspace(0, 2 * np.pi, N, endpoint=False)

w = np.sin(x)
sin = np.fft.fft(np.sin(x))
e = np.cos(x)
cos = np.fft.fft(np.cos(x))

l = np.linspace(-5, 5, N, endpoint=False)

e1 = np.fft.fft(np.exp(- l ** 2))
r = np.exp(- l ** 2)
e2 = np.fft.fft(np.exp(-(5 - l) ** 2) + np.exp(-(5 + l) ** 2))
t = np.exp(-(5 - l) ** 2) + np.exp(-(5 + l) ** 2)

# Real
plt.figure()
#plt.plot(w, 'o--', markersize=2, label="bez furjē sin")
plt.plot(sin.real, label="sin")
plt.xlabel("N")
plt.ylabel("real")
plt.legend()
plt.grid()
plt.show()

# Imag
plt.figure()
#plt.plot(w, 'o--', markersize=2, label="bez furjē sin")
plt.plot(sin.imag, label="sin")
plt.xlabel("N")
plt.ylabel("imag")
plt.legend()
plt.grid()
plt.show()

pld.na 
# Real
plt.figure()
plt.plot(cos.real, label="cos")
#plt.plot(e, 'o--', markersize=2, label="bez furjē cos")
plt.xlabel("N")
plt.ylabel("real")
plt.legend()
plt.grid()
plt.show()

# Imag
plt.figure()
plt.plot(cos.imag, label="cos")
#plt.plot(e, 'o--', markersize=2, label="bez furjē cos")
plt.xlabel("N")
plt.ylabel("imag")
plt.legend()
plt.grid()
plt.show()

# Real
plt.figure()
plt.plot(e1.real, label="e1")
#plt.plot(r, 'o--', markersize=2, label="bez furjē e1")
plt.xlabel("N")
plt.ylabel("real")
plt.legend()
plt.grid()
plt.show()

# Imag
plt.figure()
plt.plot(e1.imag, label="e1")
#plt.plot(r, 'o--', markersize=2, label="bez furjē e1")
plt.xlabel("N")
plt.ylabel("imag")
plt.legend()
plt.grid()
plt.show()

# Real
plt.figure()
plt.plot(e2.real, label="e2")
#plt.plot(t, 'o--', markersize=2, label="bez furjē e2")
plt.xlabel("N")
plt.ylabel("real")
plt.legend()
plt.grid()
plt.show()

# Imag
plt.figure()
plt.plot(e2.imag, label="e2")
#plt.plot(t, 'o--', markersize=2, label="bez furjē e2")
plt.xlabel("N")
plt.ylabel("imag")
plt.legend()
plt.grid()
plt.show()
