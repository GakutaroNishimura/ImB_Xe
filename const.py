import numpy as np
import pandas as pd

# # 129Xe
# #E_0 = -50.00
# E_0 = -21.81
# E_1 = 9.570
# E_2 = 90.53
# E_He = -211*10**3

# gnGamma_0 = 23.09*10**(-3)*np.sqrt(np.abs(E_0))
# #gnGamma_0 = 1.00
# gnGamma_1 = 9.24*10**(-3)
# gnGamma_2 = 5.6*10**(-3)
# nGamma_He = 954.4

# Gammagamma_0 = 121*10**(-3)
# #Gammagamma_0 = 152*10**(-3)
# Gammagamma_1 = 116*10**(-3)
# Gammagamma_He = 0.0119

# # g = (2J+1)/(2(2I+1))
# g0 = 1/4
# g1 = 3/4

# R_0 = 4.7775*10**(-15)
# R_He = 3.0*10**(-15)

# mu_N = -9.6623651*10**(-27) #[J/T]
# mu_N = mu_N/(1.62*10**(-19)) #[eV/T]
# m_N = 9.3956542052*10**8 #[eV]
# h_bar = 6.582119569*10**(-16) #[eV s]
# c = 2.99792458*10**8 #[m/s]

# # p_129 = 3*0.264 #[atm]
# p_129 = 1.0 #[atm]
# num_129 = p_129*5.88*10**3*6.02*10**23/131.3 #[m^-3]
# # d_cell = 5*10**(-2)
# d_cell = 10*10**(-2)

# p_He = 0.7
# rho_d_He = 20*2.686*10**19*10**4 #[m^-2]
# dSigma_He = 276.4/10**28 #[m^2]

He131CSList = []
He129CSList = []
Xe131CSList = []
Xe129CSList = []
dfHe = pd.read_table("./He_totalCS.txt", names=["energy", "CS"], sep="    ", skiprows=1, engine="python")
dfXe129 = pd.read_table("./Xe129_totalCS.txt", names=["energy", "CS"], sep="    ", skiprows=1, engine="python")
dfXe131 = pd.read_table("./Xe131_totalCS.txt", names=["energy", "CS"], sep="    ", skiprows=1, engine="python")
for i in range(len(dfHe.energy)):
    if 2.7 < dfHe.energy[i] < 3.7:
        He131CSList.append(dfHe.CS[i]/10**28)
for i in range(len(dfXe131.energy)):
    if 2.7 < dfXe131.energy[i] < 3.7:
        Xe131CSList.append(dfXe131.CS[i]/10**28)

for i in range(len(dfHe.energy)):
    if 9.07 < dfHe.energy[i] < 10.07:
        He129CSList.append(dfHe.CS[i]/10**28)
for i in range(len(dfXe131.energy)):
    if 9.07 < dfXe129.energy[i] < 10.07:
        Xe129CSList.append(dfXe129.CS[i]/10**28)

He131CS = np.mean(He131CSList)
He129CS = np.mean(He129CSList)
Xe131CS = np.mean(Xe131CSList)
Xe129CS = np.mean(Xe129CSList)

# 131Xe
E_0 = -5.41
E_1 = 3.2
E_2 = 14.41
E_He = -211*10**3

gnGamma_0 = 2.93*10**(-3) #1 eVでの2gGamma.
gnGamma_1 = 3.2*10**(-7) #E_1での2gGamma.
gnGamma_2 = 2.68*10**(-1) #E_2での2gGamma.
nGamma_He = 954.4

Gammagamma_0 = 1.237*10**(-1)
#Gammagamma_0 = 152*10**(-3)
Gammagamma_1 = 1.74*10**(-1)
Gammagamma_2 = 9.35*10**(-2)
Gammagamma_He = 0.0119

# g = (2J+1)/(2(2I+1))
g0 = 3/8 #J=1
g1 = 5/8 #J=2

R_0 = 4.7808*10**(-15)
R_He = 3.0*10**(-15)

mu_N = -9.6623651*10**(-27) #[J/T]
mu_N = mu_N/(1.62*10**(-19)) #[eV/T]
m_N = 9.3956542052*10**8 #[eV]
h_bar = 6.582119569*10**(-16) #[eV s]
c = 2.99792458*10**8 #[m/s]

# p_131 = 3*0.212 #[atm]
p_131 = 300*0.84/760 #[atm] for Mike
# p_131 = 10.0
num_131 = p_131*5.88*10**3*6.02*10**23/131.3 #[m^-3]
d_cell = 5*10**(-2)
# d_cell = 25*10**(-2)
# d_cell = 15*10**(-2)
# d_cell = 2*25.4*10**(-3) # Mike

p_He = 0.7
# p_He = 1.0
# rho_d_He = 20*2.686*10**19*10**4 #[m^-2]
rho_d_He = 60*2.686*10**19*10**4 #[m^-2]
dSigma_He = 276.4/10**28 #[m^2]