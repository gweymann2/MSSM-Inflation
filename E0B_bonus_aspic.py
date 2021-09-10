CConst = mp.euler + mp.log(2) - mp.mpf('2')
Cconst = CConst
C = CConst


def norm_eps3(x, alpha, phi0):
    return (mp.mpf('60') * (mp.mpf('1') + alpha * x ** 4 * (-mp.mpf('2') + x ** mp.mpf('4'))) * (
                mp.mpf('225') + alpha * x ** mp.mpf('4') * (-mp.mpf('1350') + x ** 4 * (mp.mpf('3915') + alpha * (
                    -mp.mpf('2100') + mp.mpf('20') * (mp.mpf('81') - mp.mpf('10') * alpha) * x ** 4 + mp.mpf('15') * (
                        -mp.mpf('99') + mp.mpf('20') * alpha) * x ** mp.mpf('8') + mp.mpf('90') * alpha * x ** mp.mpf(
                '12') + mp.mpf('9') * alpha * x ** 16))))) / (x ** 2 * (
                mp.mpf('15') + alpha * x ** mp.mpf('4') * (-mp.mpf('10') + mp.mpf('3') * x ** 4)) ** 2 * (
                                                                          mp.mpf('15') + alpha * x ** 4 * (
                                                                              mp.mpf('40') + x ** 4 * (
                                                                                  -mp.mpf('78') + alpha * (
                                                                                      mp.mpf('20') + mp.mpf(
                                                                                  '3') * x ** 8))))) / phi0 ** 2


def find_reheat_rrad_anyorder(nuStar, calFplusNuEnd, epsStarVec, Vstar):
    return nuStar - calFplusNuEnd + mp.mpf('0.25') * mp.log(
        (mp.mpf('9') - mp.mpf('3') * epsStarVec[0]) / (epsStarVec[0] * Vstar)) + mp.mpf('0.25') * mp.log(
        slowroll_corrections(epsStarVec))


def slowroll_corrections(epsV):
    return hubbleflow_corrections(slowroll_to_hubble(epsV))


def slowroll_to_hubble(epsV):
    epsVnp1 = mp.mpf('0')
    neps = len(epsV)

    if neps == 1 or neps == 2:
        return epsV
    elif neps == 3:
        a = epsV[0] * (mp.mpf('1') - epsV[1] / 3)
        b = epsV[1] * (mp.mpf('1') - epsV[1] / 6 - epsV[2] / 3)
        c = epsV[2] * (mp.mpf('1') - epsV[1] / 3 - epsVnp1 / 3)
        epsH = [a, b, c]
    else:
        print('slowroll_to_hubble: order not implemented!')
        return mp.mpf('nan')
    return epsH


def hubbleflow_corrections(epsH):
    neps = len(epsH)
    if neps == 1:
        hubbleflow_corrections = mp.mpf('1')

    elif neps == 2:
        hubbleflow_corrections = mp.mpf('1') - 2 * (mp.mpf('1') + Cconst) * epsH[0] - CConst * epsH[1]

    elif neps == 3:
        hubbleflow_corrections = mp.mpf('1') - 2 * (mp.mpf('1') + Cconst) * epsH[0] - CConst * epsH[1] + (
                    mp.mpf('-3') + 2 * CConst + 2 * CConst ** 2 + mp.pi ** 2 / 2) * epsH[0] ** 2 + (
                                             mp.mpf('-6') - Cconst + CConst ** 2 + mp.mpf('7') * mp.pi ** 2 / 12) * \
                                 epsH[0] * epsH[1] + (mp.mpf('-1') + CConst ** 2 / 2 + mp.pi ** 2 / 8) * epsH[
                                     1] ** 2 + (-CConst ** 2 / 2 + mp.pi ** 2 / 24) * epsH[1] * epsH[2]
    else:
        print('inverse_slow_corrections: order not implemented!')
        return mp.mpf('nan')
    return hubbleflow_corrections


def scalar_running(eps1, eps2, eps3):
    return -2 * eps1 * eps2 - eps2 * eps3


def scalar_index_NLO(eps1, eps2, eps3):
    return 1 - 2 * eps1 - eps2 - 2 * eps1 ** 2 - (3 + 2 * C) * eps1 * eps2 - C * eps2 * eps3
