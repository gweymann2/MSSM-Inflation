from E0_aspic import aspic
import pandas as pd
import mpmath as mp

nphi0 = 50
phi0_start = mp.mpf('1e17')
phi0_end = mp.mpf('1e18')

lnRrad = 0
ns_f = 0.9779

phi0B_list, xstar_list, alpha_list, mphi_list, A6_list, lambda6_list, ns_list, r_list = aspic(lnRrad, ns_f, phi0_start, phi0_end, nphi0)

data = pd.DataFrame(data={'phi0B':phi0B_list, 'xstar':xstar_list, 'alpha':alpha_list, 'mphi':mphi_list, 'A6':A6_list, 'lambda6':lambda6_list, 'ns':ns_list, 'r':r_list})
name_file = 'E2_tree09779Inst.csv'
data.to_csv(name_file)

print('\n######################################################################################################################################\n')
print('\n######################################################################################################################################\n')

lnRrad = 0
ns_f = 0.9551

phi0B_list, xstar_list, alpha_list, mphi_list, A6_list, lambda6_list, ns_list, r_list = aspic(lnRrad, ns_f, phi0_start, phi0_end, nphi0)

data = pd.DataFrame(data={'phi0B':phi0B_list, 'xstar':xstar_list, 'alpha':alpha_list, 'mphi':mphi_list, 'A6':A6_list, 'lambda6':lambda6_list, 'ns':ns_list, 'r':r_list})
name_file = 'E3_tree09551Inst.csv'
data.to_csv(name_file)

print('\n######################################################################################################################################\n')
print('\n######################################################################################################################################\n')

lnRrad = -10
ns_f = 0.9779

phi0B_list, xstar_list, alpha_list, mphi_list, A6_list, lambda6_list, ns_list, r_list = aspic(lnRrad, ns_f, phi0_start, phi0_end, nphi0)

data = pd.DataFrame(data={'phi0B':phi0B_list, 'xstar':xstar_list, 'alpha':alpha_list, 'mphi':mphi_list, 'A6':A6_list, 'lambda6':lambda6_list, 'ns':ns_list, 'r':r_list})
name_file = 'E4_tree09779m10.csv'
data.to_csv(name_file)

print('\n######################################################################################################################################\n')
print('\n######################################################################################################################################\n')

lnRrad = -10
ns_f = 0.9551

phi0B_list, xstar_list, alpha_list, mphi_list, A6_list, lambda6_list, ns_list, r_list = aspic(lnRrad, ns_f, phi0_start, phi0_end, nphi0)

data = pd.DataFrame(data={'phi0B':phi0B_list, 'xstar':xstar_list, 'alpha':alpha_list, 'mphi':mphi_list, 'A6':A6_list, 'lambda6':lambda6_list, 'ns':ns_list, 'r':r_list})
name_file = 'E5_tree09551m10.csv'
data.to_csv(name_file)

print('\n######################################################################################################################################\n')
print('\n######################################################################################################################################\n')

lnRrad = 0
ns_f = 0.9665

phi0B_list, xstar_list, alpha_list, mphi_list, A6_list, lambda6_list, ns_list, r_list = aspic(lnRrad, ns_f, phi0_start, phi0_end, nphi0)

data = pd.DataFrame(data={'phi0B':phi0B_list, 'xstar':xstar_list, 'alpha':alpha_list, 'mphi':mphi_list, 'A6':A6_list, 'lambda6':lambda6_list, 'ns':ns_list, 'r':r_list})
name_file = 'E6_tree09665Inst.csv'
data.to_csv(name_file)