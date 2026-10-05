import matplotlib.pyplot as plt
import numpy as np

L = 10
s = L / 1.3

def f(alpha):
    return L / (2 * alpha) - np.sinh(s / (2 * alpha))

def df(alpha):
    return -L / (2 * alpha**2) + s * np.cosh(s / (2 * alpha)) / (2 * alpha**2)

a = 1
b = 0
delta = 1e-10

while abs(a - b) > delta:
    b = a
    a = b - f(b) / df(b)

alpha = b 

def y(x):
    return alpha * np.cosh(x / alpha) - alpha

# 2 attēls
h_v = y(s / 2)
h_s = y(0)
h = abs(h_s - h_v)

print("h_vidus=", h_v)
print("h_sakums=", h_s)
print("h", h)

L_S = np.linspace(5, 1, 100)
#L_S1 = np.linspace(2, 1, 100)

H = []

#H1 = []

H = y((L_S / L) / 2)

#H1 = y((L_S1 / L) / 2)

#print(H)

plt.plot(L_S, H, "o-", label="h atkarība no L / s", markersize= 3, color="pink")
#plt.plot(L_S1, H1, "o-", label="h atkarība no L / s tuvu 1", markersize= 3, color="lightblue")
plt.xlabel("L_S")
plt.ylabel("h")
plt.grid(color="grey")
plt.legend()
plt.show()


# Kapēc nav h=0 bet ir h duadz.
