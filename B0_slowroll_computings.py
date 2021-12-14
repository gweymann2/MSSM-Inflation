from A_def_potentials import *

mp.dps = 500
mp.prec = 166

lnMpinGev = mp.mpf('42.334')
Mp = mp.exp(lnMpinGev)


def DV(V, phi):
    return mp.fdiv((mp.fsub(V(phi + 1), V(phi - 1))), 2)


def D2V(V, phi):
    return (V(phi + 1) + V(phi - 1) - 2 * V(phi)) / mp.mpf('1')


def eps1_(V, phi):
    return mp.fdiv(mp.power((mp.fmul((mp.mpf(1) * Mp), mp.fdiv(DV(V, phi), V(phi)))), 2), 2)


def eps2_(V, phi):
    return mp.fmul(mp.fmul(2, mp.power(Mp, 2)),
                   mp.fsub(mp.power(mp.fdiv(DV(V, phi), V(phi)), 2), mp.fdiv(D2V(V, phi), V(phi))))


def N(V, phimin, phimax):
    mp.dps = 50
    integ = mp.quad(lambda phi: V(phi) / DV(V, phi), [phimin, phimax])
    mp.dps = 500
    #     print('.', end="", flush=True)
    return -integ / Mp ** 2


def find_phi_st2(V, phi, phi_end, Pstar, lnRrad):
    kstar, lnMpcToKappa, HubbleSquareRootOf3OmegaRad, RelatDofRatio = 0.05, 130.282, 7.5437e-63, 1
    N0 = log(kstar) - lnMpcToKappa - 0.5 * log(HubbleSquareRootOf3OmegaRad) - 0.25 * log(RelatDofRatio)
    #     lnRrad = 0
    Delta_N_star = N(V, phi, phi_end)
    return -Delta_N_star + lnRrad - N0 - 0.25 * mp.log(
        9 / eps1_(V, phi) / (3 - eps1_(V, phi_end)) * V(phi_end) / V(phi)) + 0.25 * mp.log(8 * mp.pi ** 2 * Pstar)


def ns_(V, phi):
    return mp.fsub(mp.fsub(1, 2 * eps1_(V, phi)), eps2_(V, phi))


def endinf(V, start):
    return mp.findroot(lambda phi: eps1_(V, phi) - 1, 0.95 * start, tol=1e-30)


def phi_star(V, Pstar, lnRrad, phi_single):
    # extensions dispos si pas de guess, ou si on explore les potentiels non monotones
    phi_0 = mp.findroot(lambda phi: D2V(V, phi), phi_single, tol=1e-18)
    try:
        phi_0
    except:
        print('échec trouver un phi_0')
        return nan
    phi_start = phi_0
    phi_end = endinf(V, phi_start * 0.98)
    Ntot = N(V, phi_end, phi_start - 1)

    if -Ntot > 22.4:
        #         try:
        phi_star = mp.findroot(lambda phi: find_phi_st2(V, phi, phi_end, Pstar, lnRrad),
                               x0=(phi_start - 10, phi_start + 10), maxsteps=30, verbose=False, method='ridder',
                               tol=1e-25)
        return phi_star
    #         except:
    #             print('échec trouver un phi_star (malgré assez efolds...)')
    #             return mp.mpf(nan)
    else:
        print('pas assez efolds')
        return float(nan)

def D3V(V, phi):
    return (-0.5 * V(phi - 2) + V(phi - 1) - V(phi + 1) + 0.5 * V(phi + 2)) / mp.mpf('1')




def eps3_(V, phi):
    return float(mp.mpf('2') * Mp ** mp.mpf('4') / eps2f(V, phi) * (
                D3V(V, phi) * DV(V, phi) / V(phi) ** 2 - 3 * D2V(V, phi) * DV(V, phi) / V(phi) ** 3 + 2 * (
                    DV(V, phi) / V(phi)) ** 4))

def beginf(V, start):
    return mp.findroot(lambda phi: eps1f(V, phi) - 1, 1.05 * start)

def P_star(V, phi_sta):
    Vstar = V(phi_sta)
    eps1star = eps1_(V, phi_sta)
    return Vstar/(Mp**4*24*mp.pi**2*eps1star)