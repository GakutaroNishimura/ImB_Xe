import numpy as np
import matplotlib.pyplot as plt
import time
import ImRe_func as func


B_0 = 1.5*10**(-3) #[T]
E_0 = 9.57 #[eV]
Pol129 = 0.02
x = np.linspace(10*10**(-3), 10.0, 10000)
iB = np.linspace(-0.1, 0.1, 10000)
#x = np.linspace(8.2, 11.0, 10000)

A_s = []

for ix in x:
    iA_s = func.Tasymm_analyzer(ix, Pol129, B_0)
    A_s.append(iA_s)
fig, ax = plt.subplots()
#ax.plot(x, func.Omega_pseud(x), label="$\omega_p$ [Hz]")
ax.plot(x, A_s)
#ax.plot(x, func.B_pseudo(x), label="$B_p$ [T]")
#ax.plot(iB, func.Omega_prime(E_0, iB, Pol129=0.01), label="$\omega_p'$ [Hz]")
ax.set_xlabel("E [eV]")
#ax.set_xlabel("$B_0$ [T]")
#ax.set_ylabel("$\omega_p$ [Hz]")
ax.set_ylabel("$A_s$")
#ax.set_ylabel("$B_p$ [T]")
#ax.set_ylabel("$\omega_p'$ [Hz]")
plt.grid()
plt.legend()
fig.tight_layout()
plt.show()
