import numpy as np

# 129Xe
#E_0 = -50.00
E_0 = -21.81
E_1 = 9.570
E_2 = 90.53
E_He = -211*10**3

gnGamma_0 = 23.09*10**(-3)*np.sqrt(np.abs(E_0))
#gnGamma_0 = 1.00
gnGamma_1 = 9.24*10**(-3)
gnGamma_2 = 5.6*10**(-3)
nGamma_He = 954.4

Gammagamma_0 = 121*10**(-3)
#Gammagamma_0 = 152*10**(-3)
Gammagamma_1 = 116*10**(-3)
Gammagamma_He = 0.0119

# g = (2J+1)/(2(2I+1))
g0 = 1/4
g1 = 3/4

R_0 = 4.7775*10**(-15)
R_He = 3.0*10**(-15)

mu_N = -9.6623651*10**(-27) #[J/T]
mu_N = mu_N/(1.62*10**(-19)) #[eV/T]
m_N = 9.3956542052*10**8 #[eV]
h_bar = 6.582119569*10**(-16) #[eV s]
c = 2.99792458*10**8 #[m/s]

# p_129 = 3*0.264 #[atm]
p_129 = 1.0 #[atm]
num_129 = p_129*5.88*10**3*6.02*10**23/131.3 #[m^-3]
# d_cell = 5*10**(-2)
d_cell = 10*10**(-2)

p_He = 0.7
rho_d_He = 20*2.686*10**19*10**4 #[m^-2]
dSigma_He = 276.4/10**28 #[m^2]


# # 131Xe
# E_0 = -5.41
# E_1 = 3.2
# E_2 = 14.41
# E_He = -211*10**3

# gnGamma_0 = 2.93*10**(-3) #1 eVでの2gGamma.
# gnGamma_1 = 3.2*10**(-4) #E_1での2gGamma.
# gnGamma_2 = 2.68*10**(-1) #E_2での2gGamma.
# nGamma_He = 954.4

# Gammagamma_0 = 1.237*10**(-1)
# #Gammagamma_0 = 152*10**(-3)
# Gammagamma_1 = 1.74*10**(-1)
# Gammagamma_2 = 9.35*10**(-2)
# Gammagamma_He = 0.0119

# # g = (2J+1)/(2(2I+1))
# g0 = 3/8 #J=1
# g1 = 5/8 #J=2

# R_0 = 4.7808*10**(-15)
# R_He = 3.0*10**(-15)

# mu_N = -9.6623651*10**(-27) #[J/T]
# mu_N = mu_N/(1.62*10**(-19)) #[eV/T]
# m_N = 9.3956542052*10**8 #[eV]
# h_bar = 6.582119569*10**(-16) #[eV s]
# c = 2.99792458*10**8 #[m/s]

# p_131 = 3*0.212 #[atm]
# num_131 = p_131*5.88*10**3*6.02*10**23/131.3 #[m^-3]
# d_cell = 5*10**(-2)
# # d_cell = 10*10**(-2)

# p_He = 0.7
# rho_d_He = 20*2.686*10**19*10**4 #[m^-2]
# dSigma_He = 276.4/10**28 #[m^2]