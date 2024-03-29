import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import time
import ImRe_func as func
import const


B_0 = 0.4*10**(-3) #[T]
# B_0 = 0.0
B_list = [0.0, 1.5e-3]
E_0 = 9.57 #[eV]
Pol129 = 0.06
# Pol131 = 0.04*10**(-2)
Pol131 = 7.6*10**(-2) #Mike
# Pol131 = 1.0
#x = np.linspace(10*10**(-3), 10.0, 10000)
# x = np.linspace(2.0, 20.0, 10000)
# x = np.linspace(5.0, 14.0, 10000)
x = np.linspace(1.2, 5.2, 10000)
d = np.linspace(0.0, 30.0, 10000)
iB = np.linspace(-0.1, 0.1, 10000)
# x = np.linspace(5.0, 14.0, 10000)

A_s_enrich = []
A_s_enrich1 = []
A_s_enrich2 = []
A_s_nat = []
A_s_nat2 = []
A_s_nat3 = []
A_s_nat4 = []
A_s_Npol100 = []

A_s_Mike0 = []
A_s_Mike1 = []
A_s_Mike2 = []
A_s_Mike3 = []

Nlist = []


# for ix in x:
#     iA_s = func.Tasymm_analyzer(ix, Pol129, B_0, 1.0, 0.1)
#     A_s_enrich.append(iA_s)

# for ix in x:
#     iA_s = func.Tasymm_analyzer(ix, 0.02, B_0, 1.0, 0.1)
#     A_s_enrich1.append(iA_s)

# for ix in x:
#     iA_s = func.Tasymm_analyzer(ix, Pol129, B_0, 3.0*0.264, 0.05)
#     A_s_nat.append(iA_s)

# for ix in x:
#     iA_s = func.Tasymm_analyzer2(ix, Pol129, B_0, 1.0, const.d_cell)
#     A_s_enrich2.append(iA_s)

# for ix in x:
#     iA_s = func.Tasymm_analyzer2(ix, Pol131, B_0, 3.0*0.212, const.d_cell)
#     A_s_nat2.append(iA_s)
    
# for ix in x:
#     iA_s = func.Tasymm_analyzer3(ix, Pol131, B_0, 3.0*0.212, const.d_cell)
#     A_s_nat3.append(iA_s)
    
# for ix in x:
#     iA_s = func.Tasymm_analyzer3(ix, Pol131, 0.0, 3.0*0.212, const.d_cell)
#     A_s_nat4.append(iA_s)

# for ix in x:
#     iA_s = func.Tasymm_analyzer2_1(ix, Pol131, B_0, 3.0*0.212, const.d_cell)
#     A_s_nat2.append(iA_s)

# for ix in x:
#     iA_s = func.Tasymm_analyzer2(ix, Pol131, B_0, 3.0*0.212, const.d_cell)
#     jA_s = func.Tasymm_analyzer2_1(ix, Pol131, B_0, 3.0*0.212, const.d_cell)
#     A_s = iA_s - jA_s
#     A_s_nat2.append(A_s)

# for ix in x:
#     iA_s = func.Tasymm_analyzer3(ix, Pol131, B_0, 3.0*0.212, const.d_cell)
#     jA_s = func.Tasymm_analyzer3_1(ix, Pol131, B_0, 3.0*0.212, const.d_cell)
#     A_s = iA_s - jA_s
#     A_s_nat2.append(A_s)
    
for ix in x:
    # iA_s = func.Tasymm_analyzer2(ix, Pol131, B_0, 300*0.84/760, const.d_cell*0)
    # jA_s = func.Tasymm_analyzer2_1(ix, Pol131, B_0, 300*0.84/760, const.d_cell*0)
    iA_s = func.Tasymm_analyzer3(ix, Pol131, B_0, 300*0.84/760, const.d_cell*0)
    jA_s = func.Tasymm_analyzer3_1(ix, Pol131, B_0, 300*0.84/760, const.d_cell*0)
    A_s = iA_s - jA_s
    A_s_Mike0.append(A_s)

for ix in x:
    # iA_s = func.Tasymm_analyzer2(ix, Pol131, B_0, 300*0.84/760, const.d_cell)
    # jA_s = func.Tasymm_analyzer2_1(ix, Pol131, B_0, 300*0.84/760, const.d_cell)
    iA_s = func.Tasymm_analyzer3(ix, Pol131, B_0, 300*0.84/760, const.d_cell)
    jA_s = func.Tasymm_analyzer3_1(ix, Pol131, B_0, 300*0.84/760, const.d_cell)
    A_s = iA_s - jA_s
    A_s_Mike1.append(A_s)
    
for ix in x:
    # iA_s = func.Tasymm_analyzer2(ix, Pol131, B_0, 300*0.84/760, const.d_cell*2)
    # jA_s = func.Tasymm_analyzer2_1(ix, Pol131, B_0, 300*0.84/760, const.d_cell*2)
    iA_s = func.Tasymm_analyzer3(ix, Pol131, B_0, 300*0.84/760, const.d_cell*2)
    jA_s = func.Tasymm_analyzer3_1(ix, Pol131, B_0, 300*0.84/760, const.d_cell*2)
    A_s = iA_s - jA_s
    A_s_Mike2.append(A_s)
    
for ix in x:
    # iA_s = func.Tasymm_analyzer2(ix, Pol131, B_0, 300*0.84/760, const.d_cell*3)
    # jA_s = func.Tasymm_analyzer2_1(ix, Pol131, B_0, 300*0.84/760, const.d_cell*3)
    iA_s = func.Tasymm_analyzer3(ix, Pol131, B_0, 300*0.84/760, const.d_cell*3)
    jA_s = func.Tasymm_analyzer3_1(ix, Pol131, B_0, 300*0.84/760, const.d_cell*3)
    A_s = iA_s - jA_s
    A_s_Mike3.append(A_s)
    
# for ix in x:
#     iA_s = func.Tasymm_analyzer2(ix, Pol131, B_0, 3.0, const.d_cell)
#     A_s_nat2.append(iA_s)

# A_s_rate = [A_s_enrich[i]/A_s_enrich1[i] for i in range(len(A_s_enrich))]
    
# for ix in x:
#     iA_s = func.N100Pol(ix, Pol131, 10.0)
#     A_s_Npol100.append(iA_s)


# for id in d:
#     iN = func.NafterAnalyzer(3.2, Pol131, B_0, 300*0.84/760, id*1e-2)
#     # iN = iN*7*24*3600
#     # iN = 1/np.sqrt(iN)
#     Nlist.append(iN)


fig, ax = plt.subplots()

# ax.get_yaxis().get_major_formatter().set_useOffset(False)
# fig2, ax2 = plt.subplots()

# i = 0
# AsList = [[],[],[]]
# for iB0 in B_list:
#     for ix in x:
#         iA_s = func.Tasymm_analyzer2(ix, Pol129, iB0, 1.0, 0.1)
#         AsList[i].append(iA_s)
#         # print(AsList)
#     i += 1
# print(AsList)
# ax.plot(x, AsList[0], label="$B_0=%f$ [mT]"%0.0)
# ax.plot(x, AsList[1], label="$B_0=%f$ [mT]"%1.5)
# ax.plot(x, AsList[2], label="$B_0=%f$ [Hz]"%1.5)
# ax.plot(x, func.Omega_zero(B_0) + x*0, label="$\omega_0$")
# ax.plot(x, func.Omega_pseud(x, 1.0), label="$\omega_p$ [Hz]")
# ax.plot(x, func.Omega_prime(x, B_0, Pol129=0.06), label="$\omega_p'$ [Hz]")
# ax.plot(x, func.Omega_prime131(x, B_0, Pol131), label="$\omega_p'$")
# ax.plot(x, func.Omega_prime131(x, B_0, Pol131), label="$\omega_p'$")
# ax.plot(x, func.Omega_prime131_1(x, B_0, Pol131), label="$\omega_p'$ only p-wave")
# ax.plot(x, func.Omega_pseud131(x, 3.0*0.212), label="$\omega_p'$ [Hz]")
# ax.plot(x, func.Omega_pseud131_1(x, 3.0*0.212), label="$\omega_p'$ [Hz]")
#ax.plot(iB, func.Omega_prime(E_0, iB, Pol129=0.01), label="$\omega_p'$ [Hz]")
#ax.plot(x, func.B_pseudo(x), label="$B_p$ [T]")
# ax.plot(x, func.n_rotation(x, B_0, Pol129), label=r"$\theta_{\mathrm{n}}$ [deg.]")
# ax.plot(x, func.n_rotation131(x, B_0, Pol131), label=r"$\theta_{\mathrm{n}}$ [deg.]")
# ax.plot(x, func.n_rotation131(x, B_0/2, Pol131), label=r"$\theta_{\mathrm{n}}$ [deg.]")

#ax.plot(x, A_s_rate, label="rate")
# ax.plot(x, A_s_enrich, label="129 enrich : 1 atm 10 cm")
# ax.plot(x, A_s_enrich)
# ax.plot(x, A_s_enrich2)
#ax.plot(x, A_s_enrich, label="129 enrich : 1% pol")
#ax.plot(x, A_s_enrich1, label="129 enrich : 2% pol")
# ax.plot(x, A_s_nat, label="natural : 3 atm 5 cm")
# ax.plot(x, A_s_nat2, label="sinx=x")
# ax.plot(x, A_s_nat3, label="sinx")
# ax.plot(x, A_s_nat4, label="B_0=0")
# ax.plot(x, A_s_Npol100)
ax.plot(x, A_s_Mike0, label=r"$d_{\mathrm{Xe}}$ = 0 cm")
ax.plot(x, A_s_Mike1, label=r"$d_{\mathrm{Xe}}$ = 5 cm")
ax.plot(x, A_s_Mike2, label=r"$d_{\mathrm{Xe}}$ = 10 cm")
ax.plot(x, A_s_Mike3, label=r"$d_{\mathrm{Xe}}$ = 15 cm")
# ax.plot(d, Nlist)

ax.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))   # こっちを先に書くこと.
ax.ticklabel_format(style="sci", axis="y", scilimits=(0,0))   # 10^3（10の3乗）単位にする.
ax.set_xlabel("E [eV]")
# ax.set_xlabel("$d_{\mathrm{Xe}}$ [cm]")
# ax.set_xlabel("$B_0$ [T]")
ax.set_ylabel("$\omega$ [Hz]")
#ax.set_ylabel("$\omega_p$ [Hz]")
# ax.set_ylabel("$A_{\mathrm{asym}}$")
# ax.set_ylabel("Neutron Flux [n/s]")
#ax.set_ylabel("$B_p$ [T]")
#ax.set_ylabel("$\omega_p'$ [Hz]")
# ax.set_ylabel("neutron rotation [deg.]")
plt.grid()
plt.legend()
fig.tight_layout()
plt.show()
