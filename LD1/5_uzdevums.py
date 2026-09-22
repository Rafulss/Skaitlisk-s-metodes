import matplotlib.pyplot as plt
import numpy as np

# Dotā sākuma novirze
a0 = np.pi / 3

N = 10000
periodu_skaits = 2

# a = alpha
def f(a, a0):
    return 1 / np.sqrt(np.cos(a) - np.cos(a0))

# a no alpha0 līdz 0
a1 = np.linspace(a0, 0, N)

# Laika masīvs
t1 = np.zeros(N)

# Aprēķina laiku katrai a vērtībai
for i in range(N - 1):
    # Intervāla viduspunkts
    a_vid = (a1[i] + a1[i + 1]) / 2

    # Intervāla garums
    h = a1[i] - a1[i + 1]

    # t = integrālis d(alpha)/omega(alpha)
    t1[i + 1] = t1[i] + h * f(a_vid, a0) / np.sqrt(2)

# Ceturtdaļperiods un periods
T_q = t1[-1]
T = 4 * T_q

# 1. ceturtdaļa: alpha0 -> 0
t_1 = t1
a_1 = a1

# 2. ceturtdaļa: 0 -> -alpha0
t_2 = T_q + (T_q - t1[::-1])
a_2 = -a1[::-1]

# 3. ceturtdaļa: -alpha0 -> 0
t_3 = 2 * T_q + t1
a_3 = -a1

# 4. ceturtdaļa: 0 -> alpha0
t_4 = 3 * T_q + (T_q - t1[::-1])
a_4 = a1[::-1]

# Viens pilns periods
t_periods = np.concatenate((t_1, t_2[1:], t_3[1:], t_4[1:]))
a_periods = np.concatenate((a_1, a_2[1:], a_3[1:], a_4[1:]))

# Tukši masīvi vairākiem periodiem
t = np.zeros(len(t_periods))
a = np.zeros(len(a_periods))

t = t_periods
a = a_periods

# Pievieno nākamos periodus
for i in range(periodu_skaits - 1):
    t = np.concatenate((t, t_periods[1:] + (i + 1) * T))
    a = np.concatenate((a, a_periods[1:]))

print(f"T_tilde = {T:.6f}")

plt.figure(figsize=(9, 4.5))
plt.plot(t, a, "-")
plt.plot(t, a, "o", markersize=0.01)

plt.axhline(0, color="black", linewidth=0.8)

plt.xlabel(r"$\tilde{t}$")
plt.ylabel(r"$\alpha(\tilde{t})$")
#plt.title(r"Matemātiskā svārsta kustība")
plt.grid(True)
plt.tight_layout()
plt.show()