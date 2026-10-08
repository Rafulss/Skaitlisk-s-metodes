import matplotlib.pyplot as plt
import numpy as np

# Vada garums
L = 10

# Dotā attiecība L/s = 1.3
s = L / 1.3


def f(alpha):
    return L / (2 * alpha) - np.sinh(s / (2 * alpha))

def df(alpha):
    return (-L / (2 * alpha**2)
            + s * np.cosh(s / (2 * alpha)) / (2 * alpha**2))

a = 1
b = 0
delta = 1e-10

while abs(a - b) > delta:
    b = a
    a = b - f(b) / df(b)

alpha = a

def y(x):
    return alpha * np.cosh(x / alpha) - alpha

# Augstuma starpība pie L/s = 1.3
h_v = y(s / 2)
h_s = y(0)
h = abs(h_s - h_v)

print("L/s =", L / s)
print("alpha =", alpha)
print("h_vidus =", h_v)
print("h_sakums =", h_s)

# Attiecība L/s.
# Pie L/s = 1 vads ir taisns, tādēļ h = 0.
L_S = np.linspace(1, 5, 100)

H = []

for ratio in L_S:

    # Pie L/s = 1 h ir tieši 0
    if np.isclose(ratio, 1):
        H.append(0)
        continue

    # Attālums starp piekāršanas punktiem
    s_current = L / ratio

    # Vienādojums konkrētajai L/s attiecībai
    def f_current(alpha_current):
        return (L / (2 * alpha_current)
                - np.sinh(s_current / (2 * alpha_current)))

    # Alfa atrašanai izmanto bisekcijas metodi.
    # Kreisajā galā funkcijas vērtība ir negatīva,
    # bet labajā — pozitīva.
    alpha_left = 1e-10
    alpha_right = L

    while f_current(alpha_right) < 0:
        alpha_right = 2 * alpha_right

    # Bisekcijas iterācijas
    for i in range(100):
        alpha_middle = (alpha_left + alpha_right) / 2

        if f_current(alpha_middle) < 0:
            alpha_left = alpha_middle
        else:
            alpha_right = alpha_middle

    alpha_current = (alpha_left + alpha_right) / 2

    # y(0) = 0, tādēļ h = y(s/2)
    h_current = (
        alpha_current * np.cosh(s_current / (2 * alpha_current))
        - alpha_current
    )

    H.append(h_current)

plt.figure(figsize=(7, 4))

plt.plot( L_S, H, "o-", label=r"$h$ atkarībā no $L/s$", markersize=3, color="Magenta")

# Parāda arī doto gadījumu L/s = 1.3
plt.plot(L / s, h, "o", color="Midnightblue", markersize=7, label=r"Dotais gadījums: $L/s = 1.3$")

plt.xlabel(r"$L/s$", fontsize=12)
plt.ylabel(r"$h$", fontsize=12)
#plt.title(r"Augstuma starpība $h$ atkarībā no attiecības $L/s$")
plt.grid(color="grey")
plt.legend()
plt.tight_layout()
plt.show()