from mpmath import *
import scipy.optimize as spo
import matplotlib.pyplot as plt
import numpy as np
import time

lnMpinGev = 42.334
Mp = mp.exp(mpf(1) * lnMpinGev)
plt.rcParams['figure.figsize'] = 16, 10
mp.dps = 130


def Vtree(phi, mphi_GUT, A6_GUT, lambda6_GUT):
    return mp.fadd(mp.fsub(mp.fmul(mp.fmul(0.5, mp.power(mphi_GUT, 2)), mp.power(phi, 2)),
                           mp.fmul(mp.fdiv(mp.fmul(lambda6_GUT, A6_GUT), (mp.fmul(6, mp.power(Mp, 3)))),
                                   mp.power(phi, 6))),
                   mp.fmul(mp.fdiv(mp.power(lambda6_GUT, 2), mp.power(Mp, 6)), mp.power(phi, 10)))


def Vrge(mu, xi, mphi_GUT, A6_GUT, lambda6_GUT):
    mu_GUT = mp.mpf('3e16')
    b1, b2 = mp.mpf('11') / (8 * mp.pi ** 2), mp.mpf('1') / (8 * mp.pi ** 2)
    g1_GUT, g2_GUT = mp.sqrt(mp.pi / 10), mp.sqrt(mp.pi / 6)
    g1, g2 = g1_GUT / (mp.sqrt(1 - b1 * g1_GUT ** 2 * mp.log(mu / mu_GUT))), g2_GUT / (
        mp.sqrt(1 - b2 * g2_GUT ** 2 * mp.log(mu / mu_GUT)))
    m1, m2 = mphi_GUT * xi * (g1 / g1_GUT) ** 2, mphi_GUT * xi * (g2 / g2_GUT) ** 2
    mphi_GUT2 = mphi_GUT ** 2
    mphi2 = mphi_GUT2 + (mphi_GUT * xi) ** 2 - m2 ** 2 + mp.mpf('1') / mp.mpf('11') * ((mphi_GUT * xi) ** 2 - m1 ** 2)
    A6 = A6_GUT + 6 * (mphi_GUT * xi - m2) + mp.mpf('6') / mp.mpf('11') * (mphi_GUT * xi - m1)
    lambda6 = lambda6_GUT * (g2_GUT / g2) ** 6 * (g1_GUT / g1) ** (mp.mpf('6') / mp.mpf('11'))
    phi = mu
    return 0.5 * mphi2 * phi ** 2 - lambda6 * A6 / (6 * Mp ** 3) * phi ** 6 + lambda6 ** 2 * phi ** 10 / Mp ** 6


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


def find_phi_st(V, phi):
    return mp.mpf(24) * mp.pi ** mp.mpf(2) * mp.mpf(1) * fdiv(eps1(V, phi), V(phi)) * Pstar * power(Mp, 4)


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


def log_likelihood(ns):
    return -0.5 * ((ns - 0.9665) / (3 * 0.0038)) ** 2


def mcmc(mphi_start, A6_start, lambda6_start):
    mphi, A6, lam6 = mphi_start, A6_start, lambda6_start
    ns = ns_v2(lambda phi: Vrge(phi, 1, mphi, A6, lam6))
    A6_list, lam6_list, mphi_list, ns_list, time_step_list = [A6], [lam6], [mphi], [ns], []
    n = 1000
    count = 0
    for i in range(n):
        start = time.process_time()
        A6next = A6 + mpmathify('5e-16') * np.random.randn()
        lam6next = lam6 + mpmathify('5e-18') * np.random.randn()
        mphinext = mphi + mpmathify('5e-17') * np.random.randn()
        ns_next = ns_v2(lambda phi: Vrge(phi, 1, mphinext, A6next, lam6next))
        p = np.exp(log_likelihood(ns_next)) / np.exp(log_likelihood(ns))
        p_or_1 = min(p, 1)
        if np.random.uniform(0, 1) < p_or_1:
            print('\n' + str(ns) + ' to ' + str(ns_next) + ' of proba ' + str(p_or_1) + ' accepted.')
            A6, lam6, mphi, ns = A6next, lam6next, mphinext, ns_next
        else:
            print('\n' + str(ns) + ' to ' + str(ns_next) + ' of proba ' + str(p) + ' refused.')
        A6_list.append(A6)
        lam6_list.append(lam6)
        mphi_list.append(mphi)
        ns_list.append(ns)
        time_step = time.process_time() - start
        time_step_list.append(time_step)
        count += 1
        print("step " + str(count) + "/" + str(n) + " took " + str(time_step) + " sec.\n")

    return mphi_list, A6_list, lam6_list, ns_list, time_step_list



