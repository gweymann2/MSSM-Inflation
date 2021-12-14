from C_TREE import *

mp.dps = 500
mp.prec = 166

lnMpinGev = mp.mpf('42.334')
Mp = mp.exp(lnMpinGev)

def mphi_versus_phi_V_prime_eq_0f(phi, A6, lambda6, mphi_ini):
    return mp.findroot(lambda mphi: DV(lambda phi: Vrge(phi, 1, mphi, A6, lambda6), phi), mphi_ini, verbose=False,
                       tol=5e-18)


def mphi_and_phi0(A6, lambda6, mphi_ini, phi0_tree):
    single_point = mp.findroot(
        lambda phi: DV(lambda phi: mphi_versus_phi_V_prime_eq_0f(phi, A6, lambda6, mphi_ini), phi), phi0_tree,
        verbose=False, tol=1e-25)
    return mphi_versus_phi_V_prime_eq_0f(single_point, A6, lambda6, mphi_ini), single_point


def ns_star(mphi_start, phi0_start, A6, lambda6, Pstar, lnRrad):
    phi_sta = phi_star(lambda phi: Vrge(phi, 1, mphi_start, A6, lambda6), Pstar, lnRrad, phi0_start)
    ns_star = ns_(lambda phi: Vrge(phi, 1, mphi_start, A6, lambda6), phi_sta)
    res_P = P_star(lambda phi: Vrge(phi, 1, mphi_start, A6, lambda6), phi_sta)
    return float(ns_star)


def from_ns_find_mphi(ns, mphi_start, phi0_start, A6, lambda6, Pstar, lnRrad):
    ns_star_i = lambda i: ns_star(mphi_start * (mp.mpf('1') + i * mp.mpf('1e-24')), phi0_start, A6, lambda6, Pstar,
                                  lnRrad)
    i_tuned = mp.findroot(lambda i: ns_star_i(i) - ns, 0, verbose=False, tol=1e-12)
    mphi_tuned = mphi_start * (mp.mpf('1') + i_tuned * mp.mpf('1e-24'))
    return mphi_tuned


def normalize(mphi_nsed, A6, lambda6, phi0_start, Pstar, lnRrad):
    V_to_renormalize = lambda phi: Vrge(phi, 1, mphi_nsed, A6, lambda6)
    phi0_to_renormalize = mp.findroot(lambda phi: D2V(V_to_renormalize, phi), phi0_start, tol=5e-18)
    phi_sta_to_renormalize = phi_star(V_to_renormalize, Pstar, lnRrad, phi0_to_renormalize)
    P_star_to_renormalize = P_star(V_to_renormalize, phi_sta_to_renormalize)
    P_star_cosmo = Pstar
    normalization = P_star_cosmo / P_star_to_renormalize
    mphi, A6, lambda6 = mphi_nsed * mp.sqrt(normalization), A6 * mp.sqrt(normalization), lambda6 * mp.sqrt(
        normalization)
    ns_star = ns_(lambda phi: Vrge(phi, 1, mphi, A6, lambda6), phi_sta_to_renormalize)
    return mphi, A6, lambda6, phi0_to_renormalize, phi_sta_to_renormalize


def rge_point_generator2(A6, lambda6, mphi_tree, phi0_tree, ns_input, Pstar, lnRrad):
    # 1)
    mphi_infl, phi0_infl = mphi_and_phi0(A6, lambda6, mphi_tree, phi0_tree)

    # 2)
    alpha = aspic(lnRrad, ns_input, [phi0_tree], Pstar)[2][0]
    mphi_start = mphi_infl / mp.sqrt(alpha)
    mphi_nsed = from_ns_find_mphi(ns_input, mphi_start, phi0_infl, A6, lambda6, Pstar, lnRrad)

    # 3)
    mphi_renormalized, A6_renormalized, lambda6_renormalized, phi0_to_renormalize, phi_sta_to_renormalize = normalize(
        mphi_nsed, A6, lambda6, phi0_infl, Pstar, lnRrad)

    V_renormalized = lambda phi: Vrge(phi, 1, mphi_renormalized, A6_renormalized, lambda6_renormalized)
    phi_sta_renormalized = phi_sta_to_renormalize

    ns_out = ns_(V_renormalized, phi_sta_renormalized)
    Pstar_out = P_star(V_renormalized, phi_sta_renormalized)
    phi0 = mp.findroot(lambda phi: D2V(V_renormalized, phi), phi0_to_renormalize, tol=5e-18)
    r = 16 * eps1_(V_renormalized, phi_sta_renormalized)

    return mphi_renormalized, A6_renormalized, lambda6_renormalized, ns_out, Pstar_out, phi0, r, phi_sta_renormalized


def contour2(phi0_input_list, lnRrad, ns_input, As_input):
    compatible_mphi, compatible_A6, compatible_lambda6, compatible_ns, compatible_Ps, compatible_phi0, compatible_r, compatible_phi_star = [], [], [], [], [], [], [], []
    for i, phi0_input in enumerate(phi0_input_list):
        #         print('step '+str(i+1)+'/'+str(len(phi0_input_list))+' (scale = '+"%.5e"%phi0_input+' GeV):\n')
        phi0B_list, phistar_list, alpha_list, mphi_list, A6_list, lambda6_list, ns_list, r_list = aspic(lnRrad,
                                                                                                        ns_input,
                                                                                                        [phi0_input],
                                                                                                        As_input)
        A6, lambda6, mphi_tree, phi0_tree = A6_list[0], lambda6_list[0], mphi_list[0], phi0B_list[0]

        #         print('')

        mphi_renormalized, A6_renormalized, lambda6_renormalized, ns_out, Ps_out, phi0, r, phi_starr = rge_point_generator2(
            A6, lambda6, mphi_tree, phi0_tree, ns_input, As_input, lnRrad)

        compatible_mphi.append(mphi_renormalized)
        compatible_A6.append(A6_renormalized)
        compatible_lambda6.append(lambda6_renormalized)
        compatible_ns.append(ns_out)
        compatible_Ps.append(Ps_out)
        compatible_phi0.append(phi0)
        compatible_r.append(r)
        compatible_phi_star.append(phi_starr)
        #         print('\n-------------------------------------------------------------------------------------------------------------------------------\n')
        print(end='.')
    return compatible_mphi, compatible_A6, compatible_lambda6, compatible_ns, compatible_Ps, compatible_phi0, compatible_r, compatible_phi_star
