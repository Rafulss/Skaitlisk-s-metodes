import matplotlib.pyplot as plt
import numpy as np

T = 1.3 * 2 * np.pi
omega = 2 * np.pi / T
omega_0 = 1.3 * omega

theta = 1
N = 100

n = np.arange(-N, N + 1)
t = np.linspace(0, T, len(n), endpoint=False)

k = 10

error = 1e-18

A = np.where(np.abs(n) == 1, theta / 2, 0).astype(complex)

for i in range(k):
    alfa = np.real(np.fft.ifft(np.fft.ifftshift(A)) * len(t))

    B = np.fft.fftshift(np.fft.fft(np.sin(alfa))) / len(t)

    A_jauns = np.zeros_like(A, dtype=complex)
    mask = n != 0
    A_jauns[mask] = (omega_0**2 / (omega**2 * n[mask] ** 2)) * B[mask]

    if np.max(np.abs(A_jauns - A)) < error:
        A = A_jauns
        break
    A = A_jauns

plt.plot(t, alfa)
plt.xlabel("Laiks t")
plt.ylabel("Leņķis alfa(t)")
plt.grid()
plt.show()