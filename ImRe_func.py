import math
import numpy as np
import const

def Gamma(E):
    gnGamma_0 = const.gnGamma_0*np.sqrt(E/np.abs(const.E_0))
    gnGamma_1 = const.gnGamma_1*np.sqrt(E/const.E_1)

    nGamma_0 = gnGamma_0/(2*const.g0)
    nGamma_1 = gnGamma_1/(2*const.g1)

    Gamma_0 = nGamma_0 + const.Gammagamma_0
    Gamma_1 = nGamma_1 + const.Gammagamma_1

    return nGamma_0, nGamma_1, Gamma_0, Gamma_1


#cross section is in [barn=10^-28 m^2]
def sigma_A(E):
    nGamma_0, nGamma_1, Gamma_0, Gamma_1 = Gamma(E)
    k = 0.6947*np.sqrt(E*10**3)*10**10
    return 10**28*(-(np.pi/(2*k**2))*(nGamma_0*(2*k*(E-const.E_0)*const.R_0-Gamma_0/2)/((E-const.E_0)**2+(Gamma_0/2)**2) + 3*nGamma_1*(2*k*(E-const.E_1)*const.R_0-Gamma_1/2)/((E-const.E_1)**2+(Gamma_1/2)**2)) + 4*np.pi*const.R_0**2)
    #return 10**28*(-(np.pi/(2*k**2))*(nGamma_0*(2*k*(E-E_0)*R_0-(1-k**2*R_0**2)*Gamma_0/2)/((E-E_0)**2+(Gamma_0/2)**2) + 3*nGamma_1*(2*k*(E-E_1)*R_0-(1-k**2*R_0**2)*Gamma_1/2)/((E-E_1)**2+(Gamma_1/2)**2)) + 4*np.pi*R_0**2)
    #return 10**28*(-(np.pi/(2*k**2))*(3*nGamma_1*(2*k*(E-E_1)*R_0-Gamma_1/2)/((E-E_1)**2+(Gamma_1/2)**2)) + 3*np.pi*R_0**2)
    #return 10**28*(-(np.pi/(2*k**2))*(3*nGamma_1*(2*k*(E-E_1)*R_0-(1-k**2*R_0**2)*Gamma_1/2)/((E-E_1)**2+(Gamma_1/2)**2)) + 3*np.pi*R_0**2)

def sigma_A_delta(E):
    nGamma_0, nGamma_1, Gamma_0, Gamma_1 = Gamma(E)
    k = 0.6947*np.sqrt(E*10**3)*10**10
    return 10**28*(-(np.pi/(2*k**2))*(3*nGamma_1*(2*k*(E-const.E_1)*const.R_0-(1-k**2*const.R_0**2)*Gamma_1/2)/((E-const.E_1)**2+(Gamma_1/2)**2)) + 3*np.pi*const.R_0**2)

def sigma_B(E, P):
    nGamma_0, nGamma_1, Gamma_0, Gamma_1 = Gamma(E)
    k = 0.6947*np.sqrt(E*10**3)*10**10
    return P*10**28*(np.pi/(2*k**2))*(nGamma_0*(2*k*(E-const.E_0)*const.R_0-Gamma_0/2)/((E-const.E_0)**2+(Gamma_0/2)**2) - nGamma_1*(2*k*(E-const.E_1)*const.R_0-Gamma_1/2)/((E-const.E_1)**2+(Gamma_1/2)**2))
    #return 10**28*(np.pi/(2*k**2))*(nGamma_0*(2*k*(E-E_0)*R_0-(1-k**2*R_0**2)*Gamma_0/2)/((E-E_0)**2+(Gamma_0/2)**2) - nGamma_1*(2*k*(E-E_1)*R_0-(1-k**2*R_0**2)*Gamma_1/2)/((E-E_1)**2+(Gamma_1/2)**2))
    #return 10**28*(np.pi/(2*k**2))*(-nGamma_1*(2*k*(E-E_1)*R_0-Gamma_1/2)/((E-E_1)**2+(Gamma_1/2)**2) + np.pi*R_0**2)

def sigma_B_delta(E):
    nGamma_0, nGamma_1, Gamma_0, Gamma_1 = Gamma(E)
    k = 0.6947*np.sqrt(E*10**3)*10**10
    return 10**28*(np.pi/(2*k**2)*(-nGamma_1*(2*k*(E-const.E_1)*const.R_0-(1-k**2*const.R_0**2)*Gamma_1/2)/((E-const.E_1)**2+(Gamma_1/2)**2)) + np.pi*const.R_0**2)

def sigma_B_res(E):
    nGamma_0, nGamma_1, Gamma_0, Gamma_1 = Gamma(E)
    k = 0.6947*np.sqrt(E*10**3)*10**10
    return 10**28*(np.pi/(2*k**2)*(nGamma_1*((1-k**2*const.R_0**2)*Gamma_1/2)/((E-const.E_1)**2+(Gamma_1/2)**2)))

def sigma_B_pot(E):
    nGamma_0, nGamma_1, Gamma_0, Gamma_1 = Gamma(E)
    k = 0.6947*np.sqrt(E*10**3)*10**10
    return 10**28*(np.pi/(2*k**2)*(-nGamma_1*(2*k*(E-const.E_1)*const.R_0)/((E-const.E_1)**2+(Gamma_1/2)**2)) + np.pi*const.R_0**2)
    #return 10**28*(np.pi/(2*k**2))*(-nGamma_1*(2*k*(E-E_1)*R_0)/((E-E_1)**2+(Gamma_1/2)**2))


# Re_B is in [fm]
def Re_B(E):
    nGamma_0, nGamma_1, Gamma_0, Gamma_1 = Gamma(E)
    k = 0.6947*np.sqrt(E*10**3)*10**10
    #return 10**15*(1/(8*k)*(nGamma_0*(E-const.E_0-Gamma_0*k*const.R_0)/((E-const.E_0)**2+(Gamma_0/2)**2) - nGamma_1*(E-const.E_1-Gamma_1*k*const.R_0)/((E-const.E_1)**2+(Gamma_1/2)**2)))
    return 10**15*(1/(8*k)*(nGamma_0*(E-const.E_0+Gamma_0*k*const.R_0)/((E-const.E_0)**2+(Gamma_0/2)**2) - nGamma_1*(E-const.E_1+Gamma_1*k*const.R_0)/((E-const.E_1)**2+(Gamma_1/2)**2)))


# Omega is in [Hz]
def Omega_pseud(E):
    iRe_B = Re_B(E)
    iRe_B = iRe_B/10**15
    return ((2*np.pi*const.h_bar*const.num_129*const.c**2)/const.m_N)*iRe_B

def Omega_zero(B_0):
    return -const.mu_N*B_0/(2*const.h_bar)

def B_pseudo(E):
    iOmega = Omega_pseud(E)
    return const.h_bar*iOmega/(2*const.mu_N)

def Omega_prime(E, B_0, Pol129):
    iRe_B = Re_B(E)
    iRe_B = iRe_B/10**15
    return ((2*np.pi*const.h_bar*const.num_129*const.c**2)/const.m_N)*(Pol129*iRe_B-(const.mu_N*const.m_N*B_0)/(4*np.pi*const.h_bar**2*const.num_129*const.c**2))
    

def Transmission(E, Pol129):
    # cosh内のXe偏極率はDsigmaの中に含まれている.
    sigma_0 = sigma_A(E)/10**28
    Dsigma = sigma_B(E, Pol129)/10**28
    print(const.num_129*Dsigma*const.d_cell)
    print(math.cosh(const.num_129*Dsigma*const.d_cell))
    return math.exp(-const.num_129*sigma_0*const.d_cell)*math.cosh(const.num_129*Dsigma*const.d_cell)


def Poln_nonE():
    return -math.tanh(const.p_He*const.rho_d_He*const.dSigma_He)

def Poln(E):
    p0 = 639.2
    p1 = -71.97
    p2 = 4.636
    p3 = -0.1133
    dSigma_He = (p0 + p1*E + p2*E**2 + p3*E**3)/10**28
    #return -math.tanh(const.p_He*const.rho_d_He*const.dSigma_He)
    return -math.tanh(const.p_He*const.rho_d_He*dSigma_He)

def Tasymm(E, Pol129):
    iIm_B = sigma_B(E, Pol129)
    iIm_B = iIm_B/10**28
    #return -Poln(E)*math.tanh(iIm_B*const.num_129*const.d_cell)
    return -math.tanh(iIm_B*const.num_129*const.d_cell)


def He_sigma_A(E):
    nGamma_He = const.nGamma_He*np.sqrt(E)
    Gamma_He = nGamma_He + const.Gammagamma_He
    k = 0.6947*np.sqrt(E*10**3)*10**10
    return 10**28*(-(np.pi/(2*k**2))*(nGamma_He*(2*k*(E-const.E_He)*const.R_He-Gamma_He/2)/((E-const.E_He)**2+(Gamma_He/2)**2)) + np.pi*const.R_He**2)


def Tasymm_analyzer(E, Pol129, B_0):
    t_in_cell = const.d_cell/(437.4*math.sqrt(E*10**3))
    iIm_B = sigma_B(E, Pol129)
    iIm_B = iIm_B/10**28
    iHe_sigma = He_sigma_A(E)
    iHe_sigma = iHe_sigma/10**28
    iOmega_pseud = Pol129*Omega_pseud(E)
    iOmega_zero = Omega_zero(B_0)
    return math.tanh(const.num_129*iIm_B*const.d_cell + const.rho_d_He*const.p_He*iHe_sigma*math.sin(iOmega_pseud*t_in_cell)*math.sin(iOmega_zero*t_in_cell))*math.tanh(const.rho_d_He*const.p_He*iHe_sigma*(1 + math.cos(iOmega_pseud*t_in_cell)*math.cos(iOmega_zero*t_in_cell)))