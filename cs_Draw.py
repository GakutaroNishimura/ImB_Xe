import numpy as np
import matplotlib.pyplot as plt
import time
import ImRe_func as func

Pol129 = 1

#x = np.linspace(2.0*10**(-3), 10.0, 10000)
x = np.linspace(8.2, 11.0, 10000)

fig, ax = plt.subplots()

#ax.plot(x, func.sigma_A(x)+func.sigma_B(x, Pol129), label="A'+B'")
#ax.plot(x, func.sigma_A(x)-func.sigma_B(x, Pol129), label="A'-B'")
#ax.plot(x, (func.sigma_A(x)+func.sigma_B(x, Pol129))/(func.sigma_A(x)-func.sigma_B(x, Pol129)), label="(A'+B')/(A'-B')")
#ax.plot(x, func.sigma_A(x), label="A'")
#ax.plot(x, func.sigma_A_delta(x), label="A' delta^2")
ax.plot(x, func.sigma_B(x, Pol129), label="B'")
#ax.plot(x, func.sigma_B_delta(x), label="B' delta^2")
#ax.plot(x, func.sigma_B_res(x), label="B' res")
#ax.plot(x, func.sigma_B_pot(x), label="B' pot")


#ax.plot(x, func.Re_B(x), label="B'")
#plt.xscale('log')
#plt.yscale('log')
ax.set_xlabel("E [eV]")
ax.set_ylabel("$4\pi/k \mathrm{Im}X' [\mathrm{bn}]$")
#ax.set_ylabel("$\mathrm{Re}B' [\mathrm{fm}]$")
plt.legend()
plt.grid()
plt.show()

