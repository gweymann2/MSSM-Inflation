import numpy as np
import time
from B0_slowroll_computings import eps1f, eps2f, eps3f, nsf, rf_, phi_starf, phi_0f
from A_def_potentials import Vrge
from mpmath import *

mphi = mpmathify(
    '7446.949961062429469452587821770735257098219127860459429229360941951255418616783845167127696640252099638394991071952483716844947908')


def log_likelihood(ns):
    return -0.5 * ((ns - 0.9665) / (3 * 0.0038)) ** 2


def mcmc(n, A6, lam6):
    A6_list, lam6_list = [A6], [lam6]

    V = lambda phi: Vrge(phi, 1, mphi, A6, lam6)
    phi_0 = phi_0f(V)
    phi_star = phi_starf(V, phi_0)
    ns, r, eps1, eps2, eps3, alpha = nsf(V, phi_star), rf_(V, phi_star), eps1f(V, phi_star), eps2f(V, phi_star), eps3f(V,
                                                                                                                      phi_star), A6 ** 2 / (
                                                 40 * mphi ** 2)
    ns_list, r_list, phi_0_list, phi_star_list, eps1_list, eps2_list, eps3_list, alpha_list, time_step_list = [ns], [
        r], [phi_0], [phi_star], [eps1], [eps2], [eps3], [alpha], []

    count = 0
    for i in range(n):
        start = time.process_time()
        A6_next, lam6_next = A6 + mpmathify('1e-16') * np.random.randn(), lam6 + mpmathify('1e-18') * np.random.randn()
        V_next = lambda phi: Vrge(phi, 1, mphi, A6_next, lam6_next)

        phi_0_next = phi_0f(V_next)
        phi_star_next = phi_starf(V_next, phi_0_next)
        ns_next = nsf(V_next, phi_star_next)

        p = np.exp(log_likelihood(ns_next)) / np.exp(log_likelihood(ns))
        p_or_1 = min(p, 1)

        if np.random.uniform(0, 1) < p_or_1:
            print('\n' + str(ns) + ' to ' + str(ns_next) + ' of proba ' + str(p_or_1) + ' accepted.')
            A6, lam6, ns, V, phi_0, phi_star = A6_next, lam6_next, ns_next, V_next, phi_0_next, phi_star_next
        else:
            print('\n' + str(ns) + ' to ' + str(ns_next) + ' of proba ' + str(p) + ' refused.')

        A6_list.append(A6)
        lam6_list.append(lam6)

        r, eps1, eps2, eps3, alpha = rf_(V, phi_star), eps1f(V, phi_star), eps2f(V, phi_star), eps3f(V,
                                                                                                    phi_star), A6 ** 2 / (
                                                 40 * mphi ** 2)
        ns_list.append(ns)
        r_list.append(r)
        phi_0_list.append(phi_0)
        phi_star_list.append(phi_star)
        eps1_list.append(eps1)
        eps2_list.append(eps2)
        eps3_list.append(eps3)
        alpha_list.append(alpha)

        time_step = time.process_time() - start
        time_step_list.append(time_step)

        count += 1
        print("step " + str(count) + "/" + str(n) + " took " + str(time_step) + " sec.\n")

    return A6_list, lam6_list, ns_list, r_list, eps1_list, eps2_list, eps3_list, alpha_list, time_step_list
