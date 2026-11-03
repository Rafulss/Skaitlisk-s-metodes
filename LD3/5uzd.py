import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
# Parametri
# --------------------------------------------------
g = 9.81
l = 100

theta = 0.5
N = 10
M = 201

n = np.arange(-N, N + 1)

omega0 = np.sqrt(g / l)

# Bezdimensionālais periods un frekvence.
T = 1.3 * 2 * np.pi
omega = 2 * np.pi / T

t = np.linspace(0, T, M, endpoint=False)

An = np.where(np.abs(n) == 1, theta / 2, 0).astype(complex)

Iter = []
amplituda = []


def alfan_k(t):
    C = np.zeros(M, dtype=complex)
    C[n % M] = An

    return M * np.fft.ifft(C)


# --------------------------------------------------
# Aprēķins — sākotnējais algoritms
# --------------------------------------------------
error = 1e-8
max_iter = 10000
k = 0
error0 = np.inf

alfa_0 = alfan_k(t)

while error0 >= error and k < max_iter:
    alfa_k = alfan_k(t)

    ampl1 = np.max(np.abs(alfa_k.real))
    Iter.append(k)
    amplituda.append(ampl1)

    sin_alfa_k = np.sin(alfa_k)

    B_all = np.fft.fft(sin_alfa_k) / M
    B = B_all[n % M]

    A_new = np.zeros_like(An)
    nonzero = n != 0

    A_new[nonzero] = (
        B[nonzero] / (omega**2 * n[nonzero]**2)
    )

    error0 = np.max(np.abs(A_new - An))

    An = A_new
    k += 1

alfa_k = alfan_k(t)


# --------------------------------------------------
# Rezultātu izdrukas
# --------------------------------------------------
gala_amplituda = np.max(np.abs(alfa_k.real))

x = np.asarray(Iter)
y = np.asarray(amplituda)

izmainas = np.abs(np.diff(y))
i_straujak = int(np.argmax(izmainas))
straujak_indices = np.array([i_straujak, i_straujak + 1])

print(f"Sākuma tuvinājuma amplitūda: {y[0]:.6f} rad")
print(f"Gala amplitūda: {gala_amplituda:.10f} rad")
print(f"Veikto iterāciju skaits: {k}")
print(f"Pēdējā koeficientu izmaiņa: {error0:.3e}")
print(f"Pieļaujamā koeficientu izmaiņa: {error:.1e}")

if error0 < error:
    print("Konverģences nosacījums ir izpildīts.")
else:
    print(
        "Sasniegts iterāciju limits; "
        "konverģences nosacījums nav izpildīts."
    )

print(
    "\nVisstraujāk amplitūda mainās starp "
    f"{x[i_straujak]}. un {x[i_straujak + 1]}. iterāciju."
)
print(
    f"Amplitūda mainās no {y[i_straujak]:.6f} "
    f"līdz {y[i_straujak + 1]:.6f} rad."
)
print(f"Izmaiņas lielums: {izmainas[i_straujak]:.6f} rad")

print("\nAmplitūda katrā saglabātajā iterācijā:")
print(f"{'Iterācija':>10} {'Amplitūda (rad)':>20}")

for it, amp in zip(x, y):
    print(f"{it:>10d} {amp:>20.10f}")

# Sākotnējais cikls gala amplitūdu sarakstā neieraksta.
print(f"{k:>10d} {gala_amplituda:>20.10f}  (gala tuvinājums)")

print(
    "\nLogaritmiskajā x asī 0. iterācija ir izlaista, "
    "jo log(0) nav definēts."
)


# --------------------------------------------------
# Kopīgs grafiku noformējums
# --------------------------------------------------
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 11,
})

figure_size = (8, 4.5)

# Veselu skaitļu atzīmes lineārajai iterāciju asij.
atzimju_solis = max(1, int(np.ceil(x[-1] / 6)))
iteraciju_atzimes = np.arange(
    0, x[-1] + 1, atzimju_solis
)


# --------------------------------------------------
# 1. grafiks: sākuma tuvinājums
# --------------------------------------------------
fig1, ax1 = plt.subplots(figsize=figure_size)

ax1.plot(t, alfa_0.real, color="teal", linewidth=2)
ax1.set_xlabel(r"Bezdimensionālais laiks $\tau=\omega_0t$")
ax1.set_ylabel(r"$\alpha^{(0)}(\tau)$, rad")
ax1.set_xlim(0, T)
ax1.grid(True, alpha=0.3)

fig1.tight_layout()


# --------------------------------------------------
# 2. grafiks: gala tuvinājums
# --------------------------------------------------
fig2, ax2 = plt.subplots(figsize=figure_size)

ax2.plot(t, alfa_k.real, color="indigo", linewidth=2)
ax2.set_xlabel(r"Bezdimensionālais laiks $\tau=\omega_0t$")
ax2.set_ylabel(r"$\alpha^{(k)}(\tau)$, rad")
ax2.set_xlim(0, T)
ax2.grid(True, alpha=0.3)

fig2.tight_layout()


# --------------------------------------------------
# 3. grafiks: amplitūda atkarībā no iterācijas
# --------------------------------------------------
fig3, ax3 = plt.subplots(figsize=figure_size)

ax3.scatter(
    x,
    y,
    color="magenta",
    s=55,
    label="Amplitūda",
    zorder=3,
)

ax3.scatter(
    x[straujak_indices],
    y[straujak_indices],
    color="darkorange",
    s=75,
    label="Visstraujākā izmaiņa",
    zorder=4,
)

ax3.set_xlabel("Iterācijas numurs k")
ax3.set_ylabel("Amplitūda, rad")
ax3.set_xticks(iteraciju_atzimes)
ax3.grid(True, alpha=0.3)
ax3.legend()

fig3.tight_layout()


# --------------------------------------------------
# 4.–7. grafiks: četri atsevišķi mērogu varianti
# --------------------------------------------------
varianti = [
    ("linear", "linear"),  # 4. grafiks: lin–lin
    ("log", "linear"),     # 5. grafiks: log–lin
    ("linear", "log"),     # 6. grafiks: lin–log
    ("log", "log"),        # 7. grafiks: log–log
]

for x_scale, y_scale in varianti:
    fig, ax = plt.subplots(figsize=figure_size)

    # Logaritmiskajā x asī nevar attēlot k = 0.
    if x_scale == "log":
        mask = x > 0
    else:
        mask = np.ones(x.size, dtype=bool)

    ax.scatter(
        x[mask],
        y[mask],
        color="magenta",
        s=55,
        zorder=3,
    )

    # Izceļam attiecīgajā mērogā redzamos galapunktus.
    visible_indices = straujak_indices[mask[straujak_indices]]

    ax.scatter(
        x[visible_indices],
        y[visible_indices],
        color="darkorange",
        s=75,
        zorder=4,
    )

    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)

    x_nosaukums = "lin" if x_scale == "linear" else "log"
    y_nosaukums = "lin" if y_scale == "linear" else "log"

    ax.set_xlabel(f"Iterācijas numurs k, {x_nosaukums}")
    ax.set_ylabel(f"Amplitūda (rad), {y_nosaukums}")

    if x_scale == "linear":
        ax.set_xticks(iteraciju_atzimes)

    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()

plt.show()