import numpy as np
import matplotlib.pyplot as plt

theta = 0.5

# Šeit var mainīt harmoniku un laika punktu skaitu
N = 5
M = 26

# Bezdimensiju periods un cikliskā frekvence
T = 1.3 * 2 * np.pi
omega = 2 * np.pi / T

# Stingrāks iterāciju kritērijs, lai var pārbaudīt
# amplitūdu ar 6 zīmīgajiem cipariem
error = 1e-13
max_iter = 10000


def alfan_k(An, n, M):
    C = np.zeros(M, dtype=complex)
    C[n % M] = An

    # Reizinām ar M, jo np.fft.ifft rezultātu dala ar M
    return M * np.fft.ifft(C)


def nakamais_tuvinajums(An, n, M):
    alfa_k = alfan_k(An, n, M)

    sin_alfa_k = np.sin(alfa_k)

    B_all = np.fft.fft(sin_alfa_k) / M
    B = B_all[n % M]

    A_new = np.zeros_like(An)
    nonzero = n != 0

    A_new[nonzero] = (
        B[nonzero] / (omega**2 * n[nonzero]**2)
    )

    # n = 0 gadījumā nedrīkst dalīt ar n**2.
    # Svārstībām ap līdzsvaru izvēlamies A_0 = 0.
    A_new[n == 0] = 0

    return B, A_new


def aprekinat(N, M):
    if N < 1:
        raise ValueError("N jābūt vismaz 1.")

    if M <= 2 * N:
        raise ValueError(
            "Nepieciešams M > 2*N, lai FFT harmonikas nepārklātos."
        )

    n = np.arange(-N, N + 1)

    # Sākuma tuvinājums
    A0 = np.where(
        np.abs(n) == 1, theta / 2, 0
    ).astype(complex)

    # 4. uzdevumā prasītie B_n^(0) un A_n^(1)
    B0, A1 = nakamais_tuvinajums(A0, n, M)

    An = A0.copy()
    error0 = np.inf
    k = 0

    while error0 >= error and k < max_iter:
        B, A_new = nakamais_tuvinajums(An, n, M)

        error0 = np.max(np.abs(A_new - An))

        An = A_new
        k += 1

    if error0 >= error:
        raise RuntimeError("Nav sasniegta prasītā konverģence.")

    # Pie tau = 0 visi eksponenti ir 1.
    # Šim atrisinājumam tas ir pozitīvais galējais stāvoklis.
    amplituda = np.sum(An).real

    return n, A0, B0, A1, An, k, error0, amplituda


# Galvenais aprēķins ar izvēlētajiem N un M
n, A0, B0, A1, An, k, error0, amplituda = aprekinat(N, M)

tau = np.linspace(0, T, M, endpoint=False)

alfa_0 = alfan_k(A0, n, M)
alfa_k = alfan_k(An, n, M)


# Rezultāti manuālai ievietošanai Overleaf tabulā
print("\nPirmās iterācijas un gala koeficienti:")
print("n       B_n^(0)         A_n^(1)         A_n^(k)")

for ni in range(N + 1):
    i = np.where(n == ni)[0][0]

    print(
        f"{ni:2d}  "
        f"{B0[i].real: .10f}  "
        f"{A1[i].real: .10f}  "
        f"{An[i].real: .10f}"
    )

print("\nTabulā parādīti n >= 0.")
print("Negatīvo indeksu koeficienti sakrīt ar pozitīvajiem")
print("skaitliskās precizitātes robežās.")

print(f"\nIzmantots N = {N}, M = {M}")
print(f"Iterāciju skaits: {k}")
print(f"Pēdējā koeficientu starpība: {error0:.3e}")
print(f"Amplitūda: {amplituda:.10f} rad")
print(
    "Lielākā imaginārā daļa: "
    f"{np.max(np.abs(alfa_k.imag)):.3e}"
)


# Pārbaudām izvēlētā laika punktu skaita ietekmi
rez_2M = aprekinat(N, 2 * M)
amplituda_2M = rez_2M[-1]

print("\nLaika punktu skaita pārbaude:")
print(f"M = {M:5d}: amplitūda = {amplituda:.10f} rad")
print(f"M = {2*M:5d}: amplitūda = {amplituda_2M:.10f} rad")
print(
    "Amplitūdu starpība: "
    f"{abs(amplituda_2M - amplituda):.3e} rad"
)


# 3. uzdevuma grafiks
plt.figure(figsize=(8, 4.5))

plt.plot(
    tau, alfa_0.real, "o",
    color="teal",
    label=r"Re[$\alpha^{(0)}(\tau)$]"
)
plt.plot(
    tau, alfa_0.imag,
    color="magenta",
    label=r"Im[$\alpha^{(0)}(\tau)$]"
)

plt.xlabel(r"$\tau$")
plt.ylabel(r"$\alpha^{(0)}(\tau)$, rad")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()


# 4. uzdevuma grafiks
plt.figure(figsize=(8, 4.5))

plt.plot(
    tau, alfa_k.real,"o",
    color="teal",
    label=r"Re[$\alpha^{(k)}(\tau)$]"
)
plt.plot(
    tau, alfa_k.imag,
    color="magenta",
    label=r"Im[$\alpha^{(k)}(\tau)$]"
)

plt.xlabel(r"$\tau$")
plt.ylabel(r"$\alpha^{(k)}(\tau)$, rad")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()


# 4. uzdevuma pārbaude:
# cik harmoniku vajag amplitūdai ar 6 zīmīgajiem cipariem?
#
# Šis punktu skaits ir atsevišķs no grafiku punktu skaita.
# To arī drīkst mainīt.
M_parbaude = 4096
N_max = 15
N_ref = 40

if M_parbaude <= 2 * max(N_max, N_ref):
    raise ValueError(
        "Precizitātes pārbaudei vajag "
        "M_parbaude > 2*max(N_max, N_ref)."
    )

amplituda_ref = aprekinat(N_ref, M_parbaude)[-1]
amplituda_ref_2M = aprekinat(N_ref, 2 * M_parbaude)[-1]
amplituda_ref_30 = aprekinat(30, M_parbaude)[-1]

print("\nSalīdzinājuma rezultāta stabilitāte:")
print(
    "Starpība, palielinot N no 30 līdz 40: "
    f"{abs(amplituda_ref - amplituda_ref_30):.3e} rad"
)
print(
    "Starpība, divkāršojot M_parbaude: "
    f"{abs(amplituda_ref - amplituda_ref_2M):.3e} rad"
)

if (
    abs(amplituda_ref - amplituda_ref_30) > 1e-10
    or abs(amplituda_ref - amplituda_ref_2M) > 1e-10
):
    raise RuntimeError(
        "Salīdzinājuma rezultāts vēl nav pietiekami stabils."
    )

# Puse no pēdējā cipara vienības pie 6 zīmīgajiem cipariem
pielaujama_kluda = (
    0.5 * 10**(
        np.floor(np.log10(abs(amplituda_ref))) - 5
    )
)

N_pietiek = None

print("\nAmplitūdas atkarība no harmoniku skaita:")
print("N       amplitūda, rad       starpība pret N=40, rad")

for N_test in range(1, N_max + 1):
    amplituda_test = aprekinat(N_test, M_parbaude)[-1]
    starpiba = abs(amplituda_test - amplituda_ref)

    print(
        f"{N_test:2d}      "
        f"{amplituda_test:.10f}          "
        f"{starpiba:.3e}"
    )

    sakrit_6_cipari = (
        f"{amplituda_test:.6g}" == f"{amplituda_ref:.6g}"
    )

    if (
        N_pietiek is None
        and starpiba < pielaujama_kluda
        and sakrit_6_cipari
    ):
        N_pietiek = N_test

print(
    "\nSalīdzinājuma amplitūda ar 6 zīmīgajiem cipariem: "
    f"{amplituda_ref:.6g} rad"
)

if N_pietiek is not None:
    print(
        "Šajā skaitliskajā pārbaudē pietiek ar "
        f"|n| <= {N_pietiek}."
    )
else:
    print("Nepieciešams palielināt N_max.")