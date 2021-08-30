from mpmath import *

lnMpinGev = 42.334
Mp = mp.exp(mpf(1) * lnMpinGev)
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

# def Vapprox