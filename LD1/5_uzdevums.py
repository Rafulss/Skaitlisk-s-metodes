import matplotlib.pyplot as plt
import numpy as np

# Fizikālie parametri
g = 9.81
l = 1.0

# Dotā sākuma novirze alpha_0
a0 = np.pi / 3

N = 10000

# a = alpha
# Tā pati funkcija kā 2. uzdevumā
def f(a, a0):
    return 1 / np.sqrt(np.cos(a) - np.cos(a0))


# ---------------------------------------------------
# Aprēķinām periodu T
# ---------------------------------------------------

# Šeit integrējam no -a0 līdz a0
q = -a0
w = a0

h = (w - q) / (N - 1)

a = np.linspace(q, w, N)

# Neizmanto galapunktus, jo a = +-a0 ir singularitāte
y = f(a[1:-1], a0)

# Integrālis no -a0 līdz a0
I = h * (3/2 * y[0] + np.sum(y[1:-1]) + 3/2 * y[-1])

# No formulas:
# T = 2 * integral(-a0 līdz a0) d(alpha) / omega(alpha)
#
# omega(alpha) = sqrt(2*g/l) * sqrt(cos(alpha) - cos(a0))
#
# T ir periods sekundēs
T = 2 * np.sqrt(l / (2 * g)) * I

# Bezdimensionālais periods
T_t = T / (2 * np.pi) * np.sqrt(g / l)

print(f"alpha0 = {a0:.6f} rad")
print(f"T = {T:.6f} s")
print(f"T_tilde = {T_t:.6f}")


# ---------------------------------------------------
# Aprēķinām vienu ceturtdaļperiodu:
# alpha0 -> 0
# ---------------------------------------------------

# Veidojam alpha vērtības no a0 līdz 0
a_q = np.linspace(a0, 0, N)

# Integrēšanai izmantojam intervālu viduspunktus,
# lai nenonāktu tieši singularitātes punktā a = a0
a_mid = (a_q[:-1] + a_q[1:]) / 2

# delta alpha ir pozitīvs garums
da = a_q[:-1] - a_q[1:]

# dt = sqrt(l/(2g)) * d(alpha) / sqrt(cos(alpha)-cos(a0))
dt = np.sqrt(l / (2 * g)) * f(a_mid, a0) * da

# Laiks no t = 0 līdz T/4
t_q = np.zeros(N)
t_q[1:] = np.cumsum(dt)

# No skaitliskā integrāļa iegūtais ceturtdaļperiods
T_q = t_q[-1]

# Lai izmantotu iepriekš aprēķināto T,
# periodu sadalām četrās vienādās daļās
T_q = T / 4

# Pielāgojam t_q, lai tas beigtos tieši pie T/4
t_q = t_q / t_q[-1] * T_q


# ---------------------------------------------------
# Izveidojam vienu pilnu periodu:
#
# alpha0 -> 0 -> -alpha0 -> 0 -> alpha0
# ---------------------------------------------------

# 1. ceturtdaļa: alpha0 -> 0
t1 = t_q
a1 = a_q

# 2. ceturtdaļa: 0 -> -alpha0
t2 = T_q + t_q[1:]
a2 = -a_q[::-1][1:]

# 3. ceturtdaļa: -alpha0 -> 0
t3 = 2 * T_q + t_q[1:]
a3 = -a_q[1:]

# 4. ceturtdaļa: 0 -> alpha0
t4 = 3 * T_q + t_q[1:]
a4 = a_q[::-1][1:]

# Viens periods
t_period = np.concatenate((t1, t2, t3, t4))
a_period = np.concatenate((a1, a2, a3, a4))


# ---------------------------------------------------
# Atkārtojam vēl vienu periodu:
# grafikā būs intervāls no 0 līdz 2T
# ---------------------------------------------------

t = np.concatenate((t_period, t_period[1:] + T))
a = np.concatenate((a_period, a_period[1:]))


# ---------------------------------------------------
# Grafiks alpha(t)
# ---------------------------------------------------

plt.figure(figsize=(9, 4.5))

plt.plot(t, a, color='blue', linewidth=1.5)

plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(T, color='gray', linestyle='--', label=r'$T$')
plt.axvline(1.5 * T, color='gray', linestyle=':', label=r'$1.5T$')

plt.xlabel(r'$t$ (s)')
plt.ylabel(r'$\alpha(t)$ (rad)')
plt.title(r'Matemātiskā svārsta kustība, $\alpha_0=\pi/3$')

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()