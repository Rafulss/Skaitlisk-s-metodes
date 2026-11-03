import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

# --------------------------------------------------
# Kopīgie parametri
# --------------------------------------------------
N = 20
M = 1024

error = 1e-13
max_iter = 10000

n = np.arange(-N, N + 1)
nonzero = n != 0

# Tālāk T apzīmē periodu bezdimensionālajā laikā tau.
# Fizikālais periods: T_fiz = T / omega0.


def alfan_k(An):
    C = np.zeros(M, dtype=complex)
    C[n % M] = An

    return (M * np.fft.ifft(C)).real


def aprekinat_B(An):
    alfa = alfan_k(An)
    B_all = np.fft.fft(np.sin(alfa)) / M

    return B_all[n % M]


# --------------------------------------------------
# Furjē metode:
# 1) dota amplitūda -> aprēķina periodu;
# 2) dots periods -> aprēķina amplitūdu.
# --------------------------------------------------
def furje(theta, T_dots=None):
    An = np.where(
        np.abs(n) == 1, theta / 2, 0
    ).astype(complex)

    if T_dots is not None:
        omega2 = (2 * np.pi / T_dots)**2

    for k in range(1, max_iter + 1):
        B = aprekinat_B(An)

        D = np.zeros_like(An)
        D[nonzero] = B[nonzero] / n[nonzero]**2

        if T_dots is None:
            # Nodrošina sum(A_n) = theta.
            omega2 = np.sum(D).real / theta

            if omega2 <= 0:
                raise RuntimeError("Iegūta nepozitīva frekvences kvadrāta vērtība.")

        A_new = D / omega2
        error0 = np.max(np.abs(A_new - An))
        An = A_new

        if error0 < error:
            T = 2 * np.pi / np.sqrt(omega2)
            amplituda = np.sum(An).real

            return An, T, amplituda, k, error0

    raise RuntimeError("Nav sasniegta konverģence.")


# --------------------------------------------------
# 1LD periods: Gausa–Ležandra kvadratūra
# --------------------------------------------------
def periods_1LD(theta, punktu_skaits=64):
    x, svari = np.polynomial.legendre.leggauss(punktu_skaits)

    # Tas pats mainīgā aizvietojums kā tavā 1LD kodā.
    # Integrālfunkcija un atvasinājums ir algebriski
    # apvienoti, lai izvairītos no galapunktu singularitātes.
    kap = np.sin(theta / 2)

    integralfunkcija = 1 / np.sqrt(
        1 - kap**2 * np.sin(np.pi * x / 2)**2
    )

    # T / T0, kur T0 = 2*pi bezdimensionālajā laikā.
    T_attieciba = 0.5 * np.sum(svari * integralfunkcija)

    return 2 * np.pi * T_attieciba


# ==================================================
# 1. PUNKTS: periods atkarībā no amplitūdas
# ==================================================

# Tie paši amplitūdas punkti, kas tavā 1LD.
theta_punkti = np.arange(1, 11) * np.pi / 20

# Blīvāks amplitūdu masīvs Furjē metodes līknei.
theta_likne = np.linspace(
    theta_punkti[0], theta_punkti[-1], 60
)

T_furje_likne = np.array([
    furje(theta)[1] for theta in theta_likne
])

T_furje_punkti = np.array([
    furje(theta)[1] for theta in theta_punkti
])

T_1LD_punkti = np.array([
    periods_1LD(theta) for theta in theta_punkti
])

# Salīdzinām periodu attiecību T/T0.
T0 = 2 * np.pi

starpibas = np.abs(
    T_furje_punkti / T0 - T_1LD_punkti / T0
)

print("\n1. PUNKTS: PERIODA SALĪDZINĀJUMS")
print(
    f"{'theta (rad)':>12} "
    f"{'T/T0 Furjē':>15} "
    f"{'T/T0 1LD':>15} "
    f"{'Starpība':>12}"
)

for theta, Tf, T1, starpiba in zip(
    theta_punkti,
    T_furje_punkti,
    T_1LD_punkti,
    starpibas,
):
    print(
        f"{theta:12.6f} "
        f"{Tf / T0:15.10f} "
        f"{T1 / T0:15.10f} "
        f"{starpiba:12.3e}"
    )

print(
    "\nLielākā T/T0 starpība: "
    f"{np.max(starpibas):.3e}"
)


# ==================================================
# 2. PUNKTS: kustības salīdzinājums pie T = 1.3*T0
# ==================================================
T_dots = 1.3 * T0

# Pozitīvs sākuma tuvinājums.
# Gala amplitūdu nosaka dotais periods, nevis šī vērtība.
An, T_3LD, theta_3LD, k, error0 = furje(
    theta=0.5,
    T_dots=T_dots,
)

# Furjē risinājums uz blīva laika režģa.
omega = 2 * np.pi / T_3LD
tau_furje = np.linspace(0, T_3LD, 1500)

alfa_furje = np.real(
    np.sum(
        An[:, None]
        * np.exp(1j * n[:, None] * omega * tau_furje),
        axis=0,
    )
)


# --------------------------------------------------
# 1LD kustība no enerģijas integrāļa
# --------------------------------------------------
def kustiba_1LD(theta, punktu_skaits=10001):
    # Aizvietojums:
    # sin(alpha/2) = sin(theta/2) * cos(u).
    #
    # Tad:
    # d(tau)/du = 1 / sqrt(1 - sin²(theta/2)*cos²(u)).
    #
    # Integrālfunkcija ir gluda arī galapunktos.
    u = np.linspace(0, np.pi / 2, punktu_skaits)
    kap = np.sin(theta / 2)

    alfa_q = 2 * np.arcsin(kap * np.cos(u))

    dtau_du = 1 / np.sqrt(
        1 - kap**2 * np.cos(u)**2
    )

    tau_q = cumulative_trapezoid(
        dtau_du, u, initial=0
    )

    T_q = tau_q[-1]

    # Viens pilns periods, izmantojot simetriju.
    tau = np.concatenate((
        tau_q,
        (2 * T_q - tau_q[::-1])[1:],
        (2 * T_q + tau_q)[1:],
        (4 * T_q - tau_q[::-1])[1:],
    ))

    alfa = np.concatenate((
        alfa_q,
        (-alfa_q[::-1])[1:],
        (-alfa_q)[1:],
        alfa_q[::-1][1:],
    ))

    return tau, alfa, 4 * T_q


# Vienāda amplitūda abām metodēm.
tau_1LD, alfa_1LD, T_1LD = kustiba_1LD(theta_3LD)

# Furjē vērtības tieši 1LD laika punktos.
# Netiek pārveidots vai pielāgots 1LD periods.
alfa_furje_1LD = np.real(
    np.sum(
        An[:, None]
        * np.exp(1j * n[:, None] * omega * tau_1LD),
        axis=0,
    )
)

kustibas_starpiba = np.max(
    np.abs(alfa_furje_1LD - alfa_1LD)
)

print("\n2. PUNKTS: KUSTĪBAS SALĪDZINĀJUMS")
print(f"Kopīgā amplitūda: {theta_3LD:.10f} rad")
print(f"3LD periods tau skalā: {T_3LD:.10f}")
print(f"1LD periods tau skalā: {T_1LD:.10f}")
print(f"Periodu starpība: {abs(T_3LD - T_1LD):.3e}")
print(f"Lielākā kustības starpība: {kustibas_starpiba:.3e} rad")
print(f"Furjē iterāciju skaits: {k}")
print(f"Pēdējā koeficientu izmaiņa: {error0:.3e}")


# --------------------------------------------------
# Grafiku noformējums
# --------------------------------------------------
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 11,
})

# 1. grafiks: periods atkarībā no amplitūdas.
fig1, ax1 = plt.subplots(figsize=(8, 4.5))

ax1.plot(
    theta_likne,
    T_furje_likne / T0,
    color="teal",
    linewidth=2,
    label="3LD: Furjē metode",
)

ax1.scatter(
    theta_punkti,
    T_1LD_punkti / T0,
    color="darkorange",
    s=45,
    zorder=3,
    label="1LD: Gausa–Ležandra kvadratūra",
)

ax1.set_xlabel(r"Amplitūda $\theta$, rad")
ax1.set_ylabel(r"$T/T_0$")
ax1.grid(True, alpha=0.3)
ax1.legend()

fig1.tight_layout()


# 2. grafiks: kustības salīdzinājums.
fig2, ax2 = plt.subplots(figsize=(8, 4.5))

ax2.plot(
    tau_furje,
    alfa_furje,
    color="teal",
    linewidth=2,
    label="3LD: Furjē metode",
)

# Retinām tikai attēlojamos punktus, nevis aprēķinu.
punktu_indeksi = np.linspace(
    0, len(tau_1LD) - 1, 61
).astype(int)

ax2.scatter(
    tau_1LD[punktu_indeksi],
    alfa_1LD[punktu_indeksi],
    color="darkorange",
    s=30,
    zorder=3,
    label=r"1LD: $\tau=\int_{\alpha}^{\theta}\frac{d\beta}{\sqrt{2(\cos\beta-\cos\theta)}}$",
)

ax2.set_xlabel(r"Bezdimensionālais laiks $\tau=\omega_0t$")
ax2.set_ylabel(r"$\alpha(\tau)$, rad")
ax2.set_xlim(0, max(T_3LD, T_1LD))
ax2.grid(True, alpha=0.3)
ax2.legend()

fig2.tight_layout()

plt.show()