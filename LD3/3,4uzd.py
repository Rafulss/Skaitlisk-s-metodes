import numpy as np
import matplotlib.pyplot as plt

g = 9.81
l = 100

theta = 0.5
N = 5
n = np.arange(-N, N + 1)

omega0 = np.sqrt(g / l)
T = 1.3 * (2 * np.pi) / 1 #nondim T = T*omega0 jeb omega0 = 1
omega = 2 * np.pi / T

#defineju sakuma noacijumu kur A=0 visiem n izņemot n=+/-1
An = np.where(np.abs(n) == 1, theta / 2, 0).astype(complex)

M = 10
t = np.linspace(0, T, M, endpoint=False)

def alfan_k(t):
    C = np.zeros(M, dtype=complex)
    C[n % M] = An

    return M * np.fft.ifft(C)

    # lai izmantotu ifft jāveic parveidojumi , lai frekvences butu seciba 0, 1, .., -1
    # Reizinām ar M, jo np.fft.ifft dala rezultātu ar M
    # burtiski norakstot - np.sum(An[:, None] * np.exp(1j * n[:, None] * omega * t), axis=0 )

error = 1e-8
max_iter = 10000
k = 0
error0 = np.inf

#sakuma alfa0
alfa_0 = alfan_k(t)

while error0 >= error and k < max_iter:

    # 2. solis: A_n^(k) -> alfa^(k)(t)
    alfa_k = alfan_k(t)

    # 3. solis: atrodam B_n^(k),
    # kur sin(alfa^(k)(t)) = summa B_n exp(i*n*omega*t)
    sin_alfa_k = np.sin(alfa_k)

    # B_n = (1 / M) * sum_j sin(alfa(t_j)) * exp(-i*n*omega*t_j)
    B_all = np.fft.fft(sin_alfa_k) / M
    B = B_all[n % M]

    # 4. solis: A_n^(k+1) = B_n^(k) / (omega^2 * n^2)
    A_new = np.zeros_like(An)
    # n = 0 gadijumu nevar rekinat, jo dala ar n^2.
    nonzero = n != 0
    A_new[nonzero] = 1**2 * B[nonzero] / (omega**2 * n[nonzero]**2)

    # 5. solis: parbaudam starpibu un iterejam
    error0 = np.max(np.abs(A_new - An))

    An = A_new
    k += 1

    #print(f"iteracija = {k}, kluda = {error:.8f}")


# Gala rezultats
alfa_k = alfan_k(t)

print("\nGala Furjē koeficienti A_n:")
for ni, Ai in zip(n, An):
    print(f"n = {ni:3d}, A_n = {Ai.real:.8f} + {Ai.imag:.8f}j")

print(f"Iterāciju skaits: {k}")
print(f"Beigu kļūda: {error0:.10f}")

plt.plot(t, alfa_0.real, "o", color="teal", label="Re[$\\alpha_0(t)$]")
plt.plot(t, alfa_0.imag, color="magenta", label="Im[$\\alpha_0(t)$]")
plt.xlabel("t")
plt.ylabel("alfa0(t)")
plt.grid()
plt.legend()
plt.show()

plt.plot(t, alfa_k.real, "o", color="teal", label="Re[$\\alpha_k(t)$]")
plt.plot(t, alfa_k.imag, color="magenta", label="Im[$\\alpha_k(t)$]")
plt.xlabel("t")
plt.ylabel("alfa_k(t)")
plt.grid()
plt.legend()
plt.show()