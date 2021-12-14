from D_RGE_phi0_ns_As import *
import pandas as pd
import time

mp.dps = 500
mp.prec = 166

lnMpinGev = mp.mpf('42.334')
Mp = mp.exp(lnMpinGev)

phi0_input_list = [mp.mpf('3e15'),mp.mpf('3e16')]
# phi0_input_list = np.array([mp.mpf(str(x)) for x in np.linspace(90000000000000, 3000000000000000, 20)])
# lnRrad_list = [mp.mpf('-10'),mp.mpf('0')]
lnRrad_list = [mp.mpf('0')]
# As_input = mp.mpf('2.2030e-9')
As_input = mp.mpf('2.10310517e-9')

#ns_input_list = [mp.mpf('0.9653')]
ns_input_list = [mp.mpf('0.9653')-mp.mpf('0.0041'),mp.mpf('0.9653'),mp.mpf('0.9653')+mp.mpf('0.0041')]
direc = 'contours/'

for lnRrad in lnRrad_list:
    for ns_input in ns_input_list:
        start_time = time.process_time()
        compatible_mphi, compatible_A6, compatible_lambda6, compatible_ns, compatible_Ps, compatible_phi0, compatible_r, compatible_phi_star = contour2(
            phi0_input_list, lnRrad, ns_input, As_input)

        data = pd.DataFrame(
            data={'phi0B': compatible_phi0, 'mphi': compatible_mphi, 'A6': compatible_A6, 'lambda6': compatible_lambda6,
                  'ns': compatible_ns, 'Ps': compatible_Ps, 'r': compatible_r, 'phi*': compatible_phi_star}, dtype=str)
        name_file = 'm' + str(int(-(lnRrad))) + '_LLe_0' + str(int(10000 * ns_input)) + '.csv'
        data.to_csv(direc + name_file)
        time_step = time.process_time() - start_time
        print(str(lnRrad) + ' ' + nstr(ns_input, 5) + ' rge ' + str(time_step) + " sec.")
        start_time = time.process_time()

        phi0B_list, phistar_list, alpha_list, mphi_list, A6_list, lambda6_list, ns_list, r_list = aspic(lnRrad,
                                                                                                        ns_input,
                                                                                                        compatible_phi0,
                                                                                                        As_input)
        A6_tree, lambda6_tree, mphi_tree, phi0_tree = A6_list[0], lambda6_list[0], mphi_list[0], phi0B_list[0]
        data = pd.DataFrame(
            data={'phi0B': phi0B_list, 'mphi': mphi_list, 'A6': A6_list, 'lambda6': lambda6_list, 'ns': ns_list,
                  'Ps': [As_input] * len(phi0_input_list), 'r': r_list, 'phi*': phistar_list}, dtype=str)
        name_file = 'm' + str(int(-(lnRrad))) + '_Tree_0' + str(int(10000 * ns_input)) + '.csv'
        data.to_csv(direc + name_file)
        time_step = time.process_time() - start_time
        print(str(lnRrad) + ' ' + nstr(ns_input, 5) + ' tree ' + str(time_step) + " sec.")


