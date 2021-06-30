from mpmath import *
import numpy as np

lnMpinGev = 42.334
Mp = mp.exp(mpf(1) * lnMpinGev)
mp.dps = 130

def DV(V, phi):
    return mp.fdiv((mp.fsub(V(phi + 1), V(phi - 1))), 2)


def D2V(V, phi):
    return (V(phi + 1) + V(phi - 1) - 2 * V(phi)) / mp.mpf('1')


def eps1(V, phi):
    return mp.fdiv(mp.power((mp.fmul((mp.mpf(1) * Mp), mp.fdiv(DV(V, phi), V(phi)))), 2), 2)


def eps2(V, phi):
    return mp.fmul(mp.fmul(2, mp.power(Mp, 2)),
                   mp.fsub(mp.power(mp.fdiv(DV(V, phi), V(phi)), 2), mp.fdiv(D2V(V, phi), V(phi))))


Pstar = mp.mpf('2.2030e-9')  # Valeur utilisée pour générer les données dans Aspic

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
        9 / eps1(V, phi) / (3 - eps1(V, phi_end)) * V(phi_end) / V(phi)) + 0.25 * mp.log(8 * mp.pi ** 2 * Pstar)


def ns(V, phi):
    return mp.fsub(mp.fsub(1, 2 * eps1(V, phi)), eps2(V, phi))


def endinf(V, start):
    return mp.findroot(lambda phi: eps1(V, phi) - 1, 0.95 * start)


def beginf(V, start):
    return mp.findroot(lambda phi: eps1(V, phi) - 1, 1.05 * start)


def ns_v2(V):
    phi_0 = mp.findroot(lambda phi: D2V(V, phi), 9.3e14)
    phi_end = endinf(V, 9.3e14)
    phi_start = phi_0

    Ntot = N(V, phi_end, phi_start - 1)
    if -Ntot > 22.4:
        try:
            phi_star = mp.findroot(lambda phi: find_phi_st2(V, phi), phi_start - 5000, maxsteps=30, tol=1e-20,
                                   verbose=False)
            return float(ns(V, phi_star))
        except:
            return float(nan)
    else:
        return float(nan)
