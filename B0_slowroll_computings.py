from mpmath import *
import numpy as np

lnMpinGev = 42.334
Mp = mp.exp(mpf(1) * lnMpinGev)
mp.dps = 130


def DV(V, phi):
    return mp.fdiv((mp.fsub(V(phi + 1), V(phi - 1))), 2)


def D2V(V, phi):
    return (V(phi + 1) + V(phi - 1) - 2 * V(phi)) / mp.mpf('1')


def D3V(V, phi):
    return (-0.5 * V(phi - 2) + V(phi - 1) - V(phi + 1) + 0.5 * V(phi + 2)) / mp.mpf('1')


def eps1f(V, phi):
    return float(mp.fdiv(mp.power((mp.fmul((mp.mpf(1) * Mp), mp.fdiv(DV(V, phi), V(phi)))), 2), 2))


def eps2f(V, phi):
    return float(mp.fmul(mp.fmul(2, mp.power(Mp, 2)),
                         mp.fsub(mp.power(mp.fdiv(DV(V, phi), V(phi)), 2), mp.fdiv(D2V(V, phi), V(phi)))))


def eps3f(V, phi):
    return float(mp.mpf('2') * Mp ** mp.mpf('4') / eps2f(V, phi) * (
                D3V(V, phi) * DV(V, phi) / V(phi) ** 2 - 3 * D2V(V, phi) * DV(V, phi) / V(phi) ** 3 + 2 * (
                    DV(V, phi) / V(phi)) ** 4))


Pstar = mp.mpf('2.2030e-9')  # Valeur utilisée pour générer les données dans Aspic


def find_phi_st(V, phi):
    return mp.mpf(24) * mp.pi ** mp.mpf(2) * mp.mpf(1) * fdiv(eps1f(V, phi), V(phi)) * Pstar * power(Mp, 4)


def N(V, phimin, phimax):
    mp.dps = 50
    integ = mp.quad(lambda phi: V(phi) / DV(V, phi), [phimin, phimax])
    mp.dps = 130
    print('.', end="", flush=True)
    return -integ / Mp ** 2


def find_phi_st2(V, phi):
    kstar, lnMpcToKappa, HubbleSquareRootOf3OmegaRad, RelatDofRatio = 0.05, 130.282, 7.5437e-63, 1
    N0 = log(kstar) - lnMpcToKappa - 0.5 * log(HubbleSquareRootOf3OmegaRad) - 0.25 * log(RelatDofRatio)
    lnRrad = 0
    phi_end = endinf(V, 9.3e14)
    Delta_N_star = N(V, phi, phi_end)
    return -Delta_N_star + lnRrad - N0 - 0.25 * mp.log(
        9 / eps1f(V, phi) / (3 - eps1f(V, phi_end)) * V(phi_end) / V(phi)) + 0.25 * mp.log(8 * mp.pi ** 2 * Pstar)


def nsf(V, phi):
    return 1 - 2 * eps1f(V, phi) - eps2f(V, phi)


#     return mp.fsub(mp.fsub(1,2*eps1f(V, phi)),eps2f(V, phi))
def rf_(V, phi):
    return 16 * eps1f(V, phi)


def endinf(V, start):
    return mp.findroot(lambda phi: eps1f(V, phi) - 1, 0.95 * start)


def beginf(V, start):
    return mp.findroot(lambda phi: eps1f(V, phi) - 1, 1.05 * start)


def phi_0f(V):
    return mp.findroot(lambda phi: D2V(V, phi), 9.3e14)


def phi_starf(V, phi_0):
    #     try :
    phi_end = endinf(V, 9.3e14)
    #     except:
    #         return float(nan)
    #     try:
    #         phi_first_eps1_min = mp.findroot(lambda phi : DV(V, phi), phi_0*0.9999999)
    #         phi_start = phi_first_eps1_min
    #         print(phi_start)
    #     except:
    phi_start = phi_0

    Ntot = N(V, phi_end, phi_start - 1)
    #     print(-Ntot)
    if -Ntot > 22.4:
        try:
            phi_star = mp.findroot(lambda phi: find_phi_st2(V, phi), phi_start - 5000, maxsteps=30, tol=1e-20,
                                   verbose=False)
            return phi_star
        except:
            return float(nan)
    else:
        return float(nan)


