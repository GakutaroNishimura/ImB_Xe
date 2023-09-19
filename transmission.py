import numpy as np
import matplotlib.pyplot as plt
import time
import ImRe_func as func

Pol129 = 0.01
#x = np.linspace(10.0*10**(-3), 20.0, 10000)
x = np.linspace(9.0, 10.2, 10000)
T_list = []
T_rate = []
T_assym_E = []
T_assym = []
Pol_n = []

for ix in x:
    #T_list.append(func.Transmission(ix, Pol129))
    T_rate.append(func.Transmission(ix, Pol129)/func.Transmission(ix, 0.0))
    #T_assym.append(func.Poln_nonE()*func.Tasymm(ix, Pol129))
    #T_assym_E.append(func.Poln(ix)*func.Tasymm(ix, Pol129))
    Pol_n.append(func.Poln(ix))
fig, ax = plt.subplots()
#ax.plot(x, T_list)
ax.plot(x, T_rate)
#ax.plot(x, T_assym, label = "no E dep.")
#ax.plot(x, T_assym_E, label = "E dep.")
#ax.plot(x, [T_assym_E[i]/T_assym[i] for i in range(len(T_assym))])
#ax.plot(x, np.abs(Pol_n))
#ax.plot(x, func.B_pseudo(x), label="$B_p$ [T]")
#ax.plot(iB, func.Omega_prime(E_0, iB), label="$\omega_p'$ [Hz]")
ax.set_xlabel("E [eV]")
#ax.set_xlabel("$B_0$ [T]")
#ax.set_ylabel("$\omega_p$ [Hz]")
#ax.set_ylabel("$B_p$ [T]")
#ax.set_ylabel("$\omega_p'$ [Hz]")
#plt.yscale('log')
plt.grid()
plt.legend()
#fig.tight_layout()
plt.show()