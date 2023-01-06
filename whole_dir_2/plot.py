import pandas as pd
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from mpmath import *

mp.dps= 150

lnMpinGev = mp.mpf('42.334')
Mp = mp.exp(lnMpinGev)

def define_plot_resolution():
    fig = plt.gcf()
    DPI = fig.get_dpi()
    fig.set_size_inches(12, 8)
    ax = plt.gca()
    for tickLabel in ax.get_xticklabels()+ax.get_yticklabels():
        tickLabel.set_fontsize(29)
    ax.yaxis.label.set_size(29)
    ax.xaxis.label.set_size(29)
    ax.yaxis.offsetText.set_fontsize(29)
    ax.xaxis.offsetText.set_fontsize(29)
    return


def sci(a):
    return("%.4e"%a)
    

##########################" TREE PLOTS "##################################


tree_10_3047_9627 = pd.read_csv('tree_10_3047_9627.csv',engine='python',dtype=str)
tree_10_3047_9665 = pd.read_csv('tree_10_3047_9665.csv',engine='python',dtype=str)
tree_10_3047_9703 = pd.read_csv('tree_10_3047_9703.csv',engine='python',dtype=str)

tree_0_3033_9665 = pd.read_csv('tree_0_3033_9665.csv',engine='python',dtype=str)
tree_0_3061_9665 = pd.read_csv('tree_0_3061_9665.csv',engine='python',dtype=str)

tree_0_3047_9627 = pd.read_csv('tree_0_3047_9627.csv',engine='python',dtype=str)
tree_0_3047_9665 = pd.read_csv('tree_0_3047_9665.csv',engine='python',dtype=str)
tree_0_3047_9703 = pd.read_csv('tree_0_3047_9703.csv',engine='python',dtype=str)
tree_0_3047_eq = pd.read_csv('tree_0_3047_10353_lin.csv',engine='python',dtype=str)

file100 = tree_10_3047_9627
file101 = tree_10_3047_9665
file102 = tree_10_3047_9703
file110 = tree_0_3047_9627
file111 = tree_0_3047_9665
file112 = tree_0_3047_9703
file_ = tree_0_3047_eq

tree_10_3047_9627_lin = pd.read_csv('tree_10_3047_9627_lin.csv',engine='python',dtype=str)
tree_10_3047_9665_lin = pd.read_csv('tree_10_3047_9665_lin.csv',engine='python',dtype=str)
tree_10_3047_9703_lin = pd.read_csv('tree_10_3047_9703_lin.csv',engine='python',dtype=str)
tree_0_3047_9627_lin = pd.read_csv('tree_0_3047_9627_lin.csv',engine='python',dtype=str)
tree_0_3047_9665_lin = pd.read_csv('tree_0_3047_9665_lin.csv',engine='python',dtype=str)
tree_0_3047_9703_lin = pd.read_csv('tree_0_3047_9703_lin.csv',engine='python',dtype=str)
tree_0_3047_eq_lin = pd.read_csv('tree_0_3047_10353_lin.csv',engine='python',dtype=str)


def pd_to_array(file, column, dtype=float):
    if dtype == float:
        return np.array([mp.mpf(x) for x in file[column]], dtype=float)
    else:
        return np.array([mp.mpf(x) for x in file[column]])

color_gut, ls_gut, lw_gut = 'lightgray', '-', 5
color_llebp1, ls_llebp1, ls_sllebp1, lw_llebp1 = 'b', '-', '--', 1.4
color_uddbp1, ls_uddbp1, ls_suddbp1, lw_uddbp1 = 'g', '-', '--', 1.4
color_llebp2, ls_llebp2, ls_sllebp2, lw_llebp2 = 'r', '-', '--', 1.4
color_uddbp2, ls_uddbp2, ls_suddbp2, lw_uddbp2 = 'c', '-', '--',  1.4
color_tree = 'orange'

plt.figure(0)
plt.plot(pd_to_array(tree_10_3047_9665_lin, 'phi0'), pd_to_array(tree_10_3047_9665_lin, 'mphi_phi0')/1000,'blue')
plt.plot(pd_to_array(tree_0_3047_9665_lin, 'phi0'), pd_to_array(tree_0_3047_9665_lin, 'mphi_phi0')/1000,'red')
plt.fill_between(pd_to_array(tree_10_3047_9627_lin, 'phi0'), pd_to_array(tree_10_3047_9627_lin, 'mphi_phi0')/1000, pd_to_array(tree_10_3047_9703_lin, 'mphi_phi0')/1000, color = 'blue', alpha=0.2, label=r'$lnR_{rad} = -10$ ; $\overline{n_s}\pm \sigma_{n_s}$')
plt.fill_between(pd_to_array(tree_0_3047_9627_lin, 'phi0'), pd_to_array(tree_0_3047_9627_lin, 'mphi_phi0')/1000, pd_to_array(tree_0_3047_9703_lin, 'mphi_phi0')/1000, color='red', alpha=0.2, label=r'$lnR_{rad} = 0$ ; $\overline{n_s}\pm \sigma_{n_s}$')
#plt.plot(pd_to_array(tree_0_3047_eq_lin, 'phi0'), pd_to_array(tree_0_3047_eq_lin, 'mphi_phi0')/1000,'black',  linestyle='--', lw=2.4, label = r'$lnR_{rad} = 0$ ; $n_{s,eq}$')
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'$m_{\phi}$ (TeV)')
plt.xlim(min(pd_to_array(tree_10_3047_9627_lin, 'phi0')), 8e15)
plt.legend(fontsize=25)
plt.ylim(1.8e-1,8e2)
plt.grid()
define_plot_resolution()

plt.locator_params(axis='y', nbins=6)
plt.locator_params(axis='x', nbins=6)
plt.savefig('save_plots/mphi_phi0_tree_reheat_lin.pdf', bbox_inches='tight')
print('mphi_phi0_tree_reheat_lin.pdf saved')

plt.figure(1)
plt.plot(pd_to_array(file101, 'phi0'), pd_to_array(file101, 'mphi_phi0')/1000,'blue')
plt.plot(pd_to_array(file111, 'phi0'), pd_to_array(file111, 'mphi_phi0')/1000,'red')
plt.fill_between(pd_to_array(file100, 'phi0'), pd_to_array(file100, 'mphi_phi0')/1000, pd_to_array(file102, 'mphi_phi0')/1000, color = 'blue', alpha=0.2, label=r'$lnR_{rad} = -10$ ; $\overline{n_s}\pm \sigma_{n_s}$')
plt.fill_between(pd_to_array(file110, 'phi0'), pd_to_array(file110, 'mphi_phi0')/1000, pd_to_array(file112, 'mphi_phi0')/1000, color='red', alpha=0.2, label=r'$lnR_{rad} = 0$ ; $\overline{n_s}\pm \sigma_{n_s}$')
#plt.plot(pd_to_array(file_, 'phi0'), pd_to_array(file_, 'mphi_phi0')/1000,'black',  linestyle='--', lw=2.4, label = r'$lnR_{rad} = 0$ ; $n_{s,eq}$')
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'$m_{\phi}$ (TeV)')
plt.xlim(1e14,3.5e16)
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.legend(fontsize=25)
plt.ylim(1.8e-1,2e4)
plt.grid()
define_plot_resolution()
plt.loglog()

plt.savefig('save_plots/mphi_phi0_tree_reheat.pdf', bbox_inches='tight')
print('mphi_phi0_tree_reheat.pdf saved')

plt.figure(2)
plt.plot(pd_to_array(file101, 'phi0'), pd_to_array(file101, 'A6_phi0')/1000,'blue')
plt.plot(pd_to_array(file111, 'phi0'), pd_to_array(file111, 'A6_phi0')/1000,'red')
plt.fill_between(pd_to_array(file100, 'phi0'), pd_to_array(file100, 'A6_phi0')/1000, pd_to_array(file102, 'A6_phi0')/1000, color = 'blue', alpha=0.2, label=r'$lnR_{rad} = -10$ ; $\overline{n_s}\pm \sigma_{n_s}$')
plt.fill_between(pd_to_array(file110, 'phi0'), pd_to_array(file110, 'A6_phi0')/1000, pd_to_array(file112, 'A6_phi0')/1000, color='red', alpha=0.2, label=r'$lnR_{rad} = 0$ ; $\overline{n_s}\pm \sigma_{n_s}$')
#plt.plot(pd_to_array(file_, 'phi0'), pd_to_array(file_, 'A6_phi0')/1000,'black',  linestyle='--',lw=2.4, label = r'$lnR_{rad} = 0$ ; $n_{s,eq}$')
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'$A_{6}$ (TeV)')
plt.legend(fontsize=25)
plt.xlim(1e14,3.5e16)
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
define_plot_resolution()
plt.loglog()
plt.ylim(1.8e-1*np.sqrt(10),2e4*np.sqrt(10))
plt.grid()

plt.savefig('save_plots/A6_phi0_tree_reheat.pdf', bbox_inches='tight')
print('A6_phi0_tree_reheat.pdf saved')

plt.figure(3)
plt.plot(pd_to_array(file101, 'phi0'), pd_to_array(file101, 'lambda6_phi0'),'blue')
plt.plot(pd_to_array(file111, 'phi0'), pd_to_array(file111, 'lambda6_phi0'),'red')
plt.fill_between(pd_to_array(file100, 'phi0'), pd_to_array(file100, 'lambda6_phi0'), pd_to_array(file102, 'lambda6_phi0'), color = 'blue', alpha=0.2, label=r'$lnR_{rad} = -10$ ; $\overline{n_s}\pm \sigma_{n_s}$')
plt.fill_between(pd_to_array(file110, 'phi0'), pd_to_array(file110, 'lambda6_phi0'), pd_to_array(file112, 'lambda6_phi0'), color='red', alpha=0.2, label=r'$lnR_{rad} = 0$ ; $\overline{n_s}\pm \sigma_{n_s}$')
#plt.plot(pd_to_array(file_, 'phi0'), pd_to_array(file_, 'lambda6_phi0'),'black',  linestyle='--',lw=2.4, label = r'$lnR_{rad} = 0$ ; $n_{s,eq}$')
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'$\lambda_{6}$')
plt.legend(fontsize=25)
plt.loglog()
plt.ylim(top= 20)
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.xlim(1e14,3.5e16)
plt.ylim(1e-5, 4)
plt.grid()
define_plot_resolution()

plt.savefig('save_plots/lambda6_phi0_tree_reheat.pdf', bbox_inches='tight')
print('lambda6_phi0_tree_reheat.pdf saved')

#plt.figure(4)
#plt.plot(pd_to_array(file101, 'mphi_phi0')/1000, pd_to_array(file101, 'A6_phi0')/1000,'blue')
#plt.plot(pd_to_array(file111, 'mphi_phi0')/1000, pd_to_array(file111, 'A6_phi0')/1000,'red')
#plt.fill_between([], [], [], color = 'blue', alpha=0.2, label=r'$ln(R_{rad}) = -10$')
#plt.fill_between([], [], [], color='red', alpha=0.2, label=r'$ln(R_{rad}) = 0$')
#plt.ylabel(r'$A_6$ (TeV)')
#plt.xlabel(r'$m_{\phi}$ (TeV)')
#plt.legend(fontsize=25)
#define_plot_resolution()
#plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
#plt.loglog()
#plt.xlim(1e14,3.5e16)
#plt.grid()

#plt.savefig('save_plots/mphi_A6_tree_reheat.pdf', bbox_inches='tight')
#print('mphi_A6_tree_reheat.pdf saved')


#plt.figure(5)
#plt.plot(pd_to_array(file101, 'phi0'), 1-pd_to_array(file101, 'alpha', 'mp'),'blue')
#plt.plot(pd_to_array(file111, 'phi0'), 1-pd_to_array(file111, 'alpha', 'mp'),'red')
#plt.fill_between(pd_to_array(file100, 'phi0'), [float(x) for x in 1-pd_to_array(file100, 'alpha', 'mp')], [float(x) for x in 1-pd_to_array(file102, 'alpha', 'mp')], color = 'blue', alpha=0.2, label=r'$lnR_{rad} = -10$ ; $\overline{n_s}\pm \sigma_{n_s}$')
#plt.fill_between(pd_to_array(file110, 'phi0'), [float(x) for x in 1-pd_to_array(file110, 'alpha', 'mp')], [float(x) for x in 1-pd_to_array(file112, 'alpha', 'mp')], color='red', alpha=0.2, label=r'$lnR_{rad} = 0$ ; $\overline{n_s}\pm \sigma_{n_s}$')
#plt.plot(pd_to_array(file_, 'phi0'), [float(x) for x in 1-pd_to_array(file_, 'alpha', 'mp')],'black',  linestyle='--',lw=2.4, label = r'$lnR_{rad} = 0$ ; $n_{s,eq}$')
#plt.loglog()
#plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
#plt.xlim(1e14,3.5e16)
#plt.ylabel(r'$\alpha$')
#plt.xlabel(r'$\phi_0$ (GeV)')
#plt.legend(fontsize=25)
#plt.grid()
#define_plot_resolution()

#plt.ylim(1e-22, 1e-12)
#plt.savefig('save_plots/alpha_phi0_tree_reheat.pdf', bbox_inches='tight')
#print('alpha_phi0_tree_reheat.pdf saved')

################################" RUNNING PLOTS "###################################

phigut = mp.mpf('3e16')
lnMpinGev = mp.mpf('42.334')
Mp = mp.exp(lnMpinGev)
pre = mp.mpf('1')/(mp.mpf('8')*mp.pi**2)
b1, b2, b3 = mp.mpf('33')/mp.mpf('5'), mp.mpf('1'), mp.mpf('-3')

mphi_tree = lambda phi, mphigut : mphigut
lambda6_tree = lambda phi, lambda6gut : lambda6gut
A6_tree = lambda phi, A6gut : A6gut

g1gut = mp.sqrt(mp.mpf('5')/mp.mpf('3'))*mp.mpf('5.45185741e-01')
g2gut = mp.mpf('6.90473022e-01')
g3gut = mp.mpf('6.84972506e-01')
m1gut = mp.mpf('1.36108022e+02')
m2gut = mp.mpf('1.14286222e+03')
m3gut = mp.mpf('8.98639714e+02')

g1 = lambda phi, g1gut : g1gut/(mp.sqrt(1-pre*b1*g1gut**2*mp.log(phi/phigut)))
g2 = lambda phi, g2gut : g2gut/(mp.sqrt(1-pre*b2*g2gut**2*mp.log(phi/phigut)))
g3 = lambda phi, g3gut : g3gut/(mp.sqrt(1-pre*b3*g3gut**2*mp.log(phi/phigut)))

M1 = lambda phi, m1gut, g1gut : m1gut*(g1(phi, g1gut)/g1gut)**mp.mpf('2')
M2 = lambda phi, m2gut, g2gut : m2gut*(g2(phi, g2gut)/g2gut)**mp.mpf('2')
M3 = lambda phi, m3gut, g3gut : m3gut*(g3(phi, g3gut)/g3gut)**mp.mpf('2')

# run from gut to phi
mphi_lle = lambda phi, mphigut, m1gut, g1gut, m2gut, g2gut : mp.sqrt(mphigut**2+(m2gut**2-M2(phi, m2gut, g2gut)**2)+mp.mpf('1')/11*(m1gut**2-M1(phi, m1gut, g1gut)**2))
A6_lle = lambda phi, A6gut, m1gut, g1gut, m2gut, g2gut  : A6gut-mp.mpf('6')*(m2gut-M2(phi, m2gut, g2gut))-mp.mpf('6')/11*(m1gut-M1(phi, m1gut, g1gut))
lambda6_lle = lambda phi, lambda6gut, g1gut, g2gut  : lambda6gut*(g2gut/g2(phi, g2gut))**mp.mpf('6')*(g1gut/g1(phi, g1gut))**(mp.mpf('6')/11)
mphi_udd = lambda phi, mphigut, m1gut, g1gut, m3gut, g3gut : mp.sqrt(mphigut**2-mp.mpf('8')/9*(m3gut**2-M3(phi, m3gut, g3gut)**2)+mp.mpf('4')/99*(m1gut**2-M1(phi, m1gut, g1gut)**2))
A6_udd = lambda phi, A6gut, m1gut, g1gut, m3gut, g3gut : A6gut+mp.mpf('16')/3*(m3gut-M3(phi, m3gut, g3gut))-mp.mpf('8')/33*(m1gut-M1(phi, m1gut, g1gut))
lambda6_udd = lambda phi, lambda6gut, g1gut, g3gut : lambda6gut*(g3gut/g3(phi, g3gut))**(mp.mpf('-16')/3)*(g1gut/g1(phi, g1gut))**(mp.mpf('8')/33)

# run from phi_start to phi_end
mphi_run_lle = lambda phi_start, phi_end, mphi_start, m1gut, g1gut, m2gut, g2gut : mp.sqrt(mphi_start**2+(M2(phi_start, m2gut, g2gut)**2-M2(phi_end, m2gut, g2gut)**2)+mp.mpf('1')/11*(M1(phi_start, m1gut, g1gut)**2-M1(phi_end, m1gut, g1gut)**2))
A6_run_lle = lambda phi_start, phi_end, A6_start, m1gut, g1gut, m2gut, g2gut : A6_start-mp.mpf('6')*(M2(phi_start, m2gut, g2gut)-M2(phi_end, m2gut, g2gut))-mp.mpf('6')/11*(M1(phi_start, m1gut, g1gut)-M1(phi_end, m1gut, g1gut))
lambda6_run_lle = lambda phi_start, phi_end, lambda6_start, g1gut, g2gut : lambda6_start*(g2(phi_start, g2gut)/g2(phi_end, g2gut))**mp.mpf('6')*(g1(phi_start, g1gut)/g1(phi_end, g1gut))**(mp.mpf('6')/11)
mphi_run_udd = lambda phi_start, phi_end, mphi_start, m1gut, g1gut, m3gut, g3gut : mp.sqrt(mphi_start**2-mp.mpf('8')/9*(M3(phi_start, m3gut, g3gut)**2-M3(phi_end, m3gut, g3gut)**2)+mp.mpf('4')/99*(M1(phi_start, m1gut, g1gut)**2-M1(phi_end, m1gut, g1gut)**2))
A6_run_udd = lambda phi_start, phi_end, A6_start, m1gut, g1gut, m3gut, g3gut : A6_start+mp.mpf('16')/3*(M3(phi_start, m3gut, g3gut)-M3(phi_end, m3gut, g3gut))-mp.mpf('8')/33*(M1(phi_start, m1gut, g1gut)-M1(phi_end, m1gut, g1gut))
lambda6_run_udd = lambda phi_start, phi_end, lambda6_start, g1gut, g3gut : lambda6_start*(g3(phi_start, g3gut)/g3(phi_end, g3gut))**(mp.mpf('6')*(-mp.mpf('8'))/9)*(g1(phi_start, g1gut)/g1(phi_end, g1gut))**(4*mp.mpf('6')/99)
mphi_run_tree = lambda phi_start, phi_end, mphi_start : mphi_start
A6_run_tree = lambda phi_start, phi_end, A6_start : A6_start
lambda6_run_tree = lambda phi_start, phi_end, lambda6_start : lambda6_start





def build_pd(name_file, infl_type = None, bp = None, name_save = None):
    if name_file[6] == "1" or bp == 1:
        g1gut = mp.sqrt(mp.mpf('5')/mp.mpf('3'))*mp.mpf('5.45185741e-01')
        g2gut = mp.mpf('6.90473022e-01')
        g3gut = mp.mpf('6.84972506e-01')
        m1gut = mp.mpf('1.36108022e+02')
        m2gut = mp.mpf('1.14286222e+03')
        m3gut = mp.mpf('8.98639714e+02')
    elif name_file[6] == '2' or bp == 2:
        g1gut = mp.mpf('0.704555660557172')
        g2gut = mp.mpf('0.690364970285155')
        g3gut = mp.mpf('0.684720653567032')
        m1gut = mp.mpf('897.774785812765')
        m2gut = mp.mpf('1789.66021839594')
        m3gut = mp.mpf('882.522487633969')
    if name_file[:3] == 'udd' or infl_type == 'udd':
        mphi_run = lambda phi_start, phi_end, mphi_start : mphi_run_udd(phi_start, phi_end, mphi_start, m1gut, g1gut, m3gut, g3gut)
        A6_run = lambda phi_start, phi_end, A6_start : A6_run_udd(phi_start, phi_end, A6_start, m1gut, g1gut, m3gut, g3gut)
        lambda6_run = lambda phi_start, phi_end, lambda6_start : lambda6_run_udd(phi_start, phi_end, lambda6_start, g1gut, g3gut)
        mphi = lambda phi, mphigut : mphi_udd(phi, mphigut, m1gut, g1gut, m3gut, g3gut)
        A6 = lambda phi, A6gut : A6_udd(phi, A6gut, m1gut, g1gut, m3gut, g3gut)
        lambda6 = lambda phi, lambda6gut : lambda6_udd(phi, lambda6gut, g1gut, g3gut)
    elif name_file[:3] == 'lle' or infl_type == 'lle':
        mphi_run = lambda phi_start, phi_end, mphi_start : mphi_run_lle(phi_start, phi_end, mphi_start, m1gut, g1gut, m2gut, g2gut)
        A6_run = lambda phi_start, phi_end, A6_start : A6_run_lle(phi_start, phi_end, A6_start, m1gut, g1gut, m2gut, g2gut)
        lambda6_run = lambda phi_start, phi_end, lambda6_start : lambda6_run_lle(phi_start, phi_end, lambda6_start, g1gut, g2gut)
        mphi = lambda phi, mphigut : mphi_lle(phi, mphigut, m1gut, g1gut, m2gut, g2gut)
        A6 = lambda phi, A6gut : A6_lle(phi, A6gut, m1gut, g1gut, m2gut, g2gut)
        lambda6 = lambda phi, lambda6gut : lambda6_lle(phi, lambda6gut, g1gut, g2gut)

    file = pd.read_csv(name_file, engine='python',dtype=str)
    file['mphi_gut'] = 0
    file['A6_gut'] = 0
    file['lambda6_gut'] = 0
    file['mphi_2'] = 0
    file['A6_2'] = 0
    file['lambda6_2'] = 0

    for i in range(len(file)):

        for key in file.keys():
            file[key].iloc[i] = mp.re(file[key].iloc[i])

        file['mphi_gut'].iloc[i] = mphi_run(file['phi0'].iloc[i], mp.mpf('3e16'), file['mphi_phi0'].iloc[i])
        file['A6_gut'].iloc[i] = A6_run(file['phi0'].iloc[i], mp.mpf('3e16'), file['A6_phi0'].iloc[i])
        file['lambda6_gut'].iloc[i] = lambda6_run(file['phi0'].iloc[i], mp.mpf('3e16'), file['lambda6_phi0'].iloc[i])
        file['mphi_2'].iloc[i] = mphi_run(file['phi0'].iloc[i], mp.mpf('2000'), file['mphi_phi0'].iloc[i])
        file['A6_2'].iloc[i] = A6_run(file['phi0'].iloc[i], mp.mpf('2000'), file['A6_phi0'].iloc[i])
        file['lambda6_2'].iloc[i] = lambda6_run(file['phi0'].iloc[i], mp.mpf('2000'), file['lambda6_phi0'].iloc[i])


    file = file.drop(columns=['Unnamed: 0'])
    file = file.reindex(columns=['phi0', 'mphi_phi0', 'A6_phi0', 'lambda6_phi0', 'mphi_2', 'A6_2', 'lambda6_2', 'mphi_gut', 'A6_gut', 'lambda6_gut', 'dv', 'phi*', 'alpha'])
    
    file['to_sort'] = 0
    for i in range(len(file)):
        file['to_sort'].iloc[i] = float(file['phi0'].iloc[i])
    file = file.sort_values(by=['to_sort'])
    file.drop(columns=['to_sort'])
    
    file['yep_or_nop'] = 'nop:'
    for i in range(len(file)):
        if mp.re(file['mphi_2'].iloc[i]) == 0 : file['yep_or_nop'].iloc[i]+=' mphi_2TeV_Im'
        if mp.re(file['mphi_gut'].iloc[i]) == 0: file['yep_or_nop'].iloc[i]+=' mphi_GUT_Im'
        if file['A6_2'].iloc[i]<0 : file['yep_or_nop'].iloc[i]+=' A6_2TeV_neg'
        if file['A6_gut'].iloc[i]<0: file['yep_or_nop'].iloc[i]+=' A6_GUT_neg'
        if file['lambda6_2'].iloc[i]>100: file['yep_or_nop'].iloc[i]+=' lambda6_2TeV_big'
        if file['lambda6_gut'].iloc[i]>100: file['yep_or_nop'].iloc[i]+=' lambda6_GUT_big'
        if file['yep_or_nop'].iloc[i] == 'nop:':file['yep_or_nop'].iloc[i] = 'yep!'

    file = file.reindex(columns=['phi0',  'yep_or_nop','mphi_phi0', 'A6_phi0', 'lambda6_phi0', 'mphi_2', 'A6_2', 'lambda6_2', 'mphi_gut', 'A6_gut', 'lambda6_gut', 'dv', 'phi*', 'alpha'])
    file = file.reset_index(drop = True)
    if name_save == None: name_save=name_file
    file.to_csv('save_file2/'+name_save)
    return file

tree_0_3047_9665_high_llebp1 = build_pd('tree_0_3047_9665.csv', infl_type = 'lle', bp = 1, name_save = 'tree_0_3047_9665_high_llebp1.csv').iloc[15:].reset_index()
tree_0_3047_9665_high_uddbp1 = build_pd('tree_0_3047_9665.csv', infl_type = 'udd', bp = 1, name_save = 'tree_0_3047_9665_high_uddbp1.csv').iloc[15:].reset_index()
tree_0_3047_9665_high_llebp2 = build_pd('tree_0_3047_9665.csv', infl_type = 'lle', bp = 2, name_save = 'tree_0_3047_9665_high_llebp2.csv').iloc[15:].reset_index()
tree_0_3047_9665_high_uddbp2 = build_pd('tree_0_3047_9665.csv', infl_type = 'udd', bp = 2, name_save = 'tree_0_3047_9665_high_uddbp2.csv').iloc[15:].reset_index()
tree_0_3047_9627_high_llebp1 = build_pd('tree_0_3047_9627.csv', infl_type = 'lle', bp = 1, name_save = 'tree_0_3047_9627_high_llebp1.csv').iloc[15:].reset_index()
tree_0_3047_9627_high_uddbp1 = build_pd('tree_0_3047_9627.csv', infl_type = 'udd', bp = 1, name_save = 'tree_0_3047_9627_high_uddbp1.csv').iloc[15:].reset_index()
tree_0_3047_9627_high_llebp2 = build_pd('tree_0_3047_9627.csv', infl_type = 'lle', bp = 2, name_save = 'tree_0_3047_9627_high_llebp2.csv').iloc[15:].reset_index()
tree_0_3047_9627_high_uddbp2 = build_pd('tree_0_3047_9627.csv', infl_type = 'udd', bp = 2, name_save = 'tree_0_3047_9627_high_uddbp2.csv').iloc[15:].reset_index()
tree_0_3047_9703_high_llebp1 = build_pd('tree_0_3047_9703.csv', infl_type = 'lle', bp = 1, name_save = 'tree_0_3047_9703_high_llebp1.csv').iloc[15:].reset_index()
tree_0_3047_9703_high_uddbp1 = build_pd('tree_0_3047_9703.csv', infl_type = 'udd', bp = 1, name_save = 'tree_0_3047_9703_high_uddbp1.csv').iloc[15:].reset_index()
tree_0_3047_9703_high_llebp2 = build_pd('tree_0_3047_9703.csv', infl_type = 'lle', bp = 2, name_save = 'tree_0_3047_9703_high_llebp2.csv').iloc[15:].reset_index()
tree_0_3047_9703_high_uddbp2 = build_pd('tree_0_3047_9703.csv', infl_type = 'udd', bp = 2, name_save = 'tree_0_3047_9703_high_uddbp2.csv').iloc[15:].reset_index()

lle_bp1_0_3047_9665 = build_pd('lle_bp1_0_3047_9665.csv')
udd_bp1_0_3047_9665 = build_pd('udd_bp1_0_3047_9665.csv')
lle_bp2_0_3047_9665 = build_pd('lle_bp2_0_3047_9665.csv')
udd_bp2_0_3047_9665 = build_pd('udd_bp2_0_3047_9665.csv')
lle_bp1_0_3047_9627 = build_pd('lle_bp1_0_3047_9627.csv')
udd_bp1_0_3047_9627 = build_pd('udd_bp1_0_3047_9627.csv')
lle_bp2_0_3047_9627 = build_pd('lle_bp2_0_3047_9627.csv')
udd_bp2_0_3047_9627 = build_pd('udd_bp2_0_3047_9627.csv')
lle_bp1_0_3047_9703 = build_pd('lle_bp1_0_3047_9703.csv')
udd_bp1_0_3047_9703 = build_pd('udd_bp1_0_3047_9703.csv')
lle_bp2_0_3047_9703 = build_pd('lle_bp2_0_3047_9703.csv')
udd_bp2_0_3047_9703 = build_pd('udd_bp2_0_3047_9703.csv')


cond1 = lle_bp1_0_3047_9665['yep_or_nop']=='yep!'
cond2 = udd_bp1_0_3047_9665['yep_or_nop']=='yep!'
cond3 = lle_bp2_0_3047_9665['yep_or_nop']=='yep!'
cond4 = udd_bp2_0_3047_9665['yep_or_nop']=='yep!'

####

def alpha_old(pd):
    pd['alpha_old'] = 0
    for i in range(len(pd)):
        pd['alpha_old'].iloc[i] = 1-mp.mpf(pd['A6_phi0'].iloc[i])**2/(20*mp.mpf(pd['mphi_phi0'].iloc[i])**2)
    return np.array([float(x) for x in pd['alpha_old']])
 
f, a = plt.subplots()
a.fill_between(np.array(tree_0_3047_9665['phi0'], dtype=float),
                 (alpha_old(tree_0_3047_9627)), 
                 (alpha_old(tree_0_3047_9703))
                 ,color=color_tree, alpha=0.3)
p1 = a.plot(np.array(tree_0_3047_9665['phi0'], dtype=float),
                 (alpha_old(tree_0_3047_9665)),color = color_tree, lw=lw_uddbp2)
p2 = a.fill(np.NaN, np.NaN, color_tree, alpha=0.3)
a.fill_between(np.array(lle_bp1_0_3047_9665[cond1]['phi0'], dtype=float),
                 (alpha_old(lle_bp1_0_3047_9627[cond1])), 
                 (alpha_old(lle_bp1_0_3047_9703[cond1])),
                 color=color_llebp1, alpha=0.3)
p3 = a.plot(np.array(lle_bp1_0_3047_9665[cond1]['phi0'], dtype=float),
                 (alpha_old(lle_bp1_0_3047_9665[cond1])),color=color_llebp1, ls = ':', lw=lw_llebp1)
p4 = a.fill(np.NaN, np.NaN, color=color_llebp1, alpha=0.3)

a.fill_between(np.array(udd_bp1_0_3047_9665[cond2]['phi0'], dtype=float),
                 (alpha_old(udd_bp1_0_3047_9627[cond2])), 
                 (alpha_old(udd_bp1_0_3047_9703[cond2])),
                 color=color_uddbp1, alpha=0.3)
p5 = a.plot(np.array(udd_bp1_0_3047_9665[cond2]['phi0'], dtype=float),
                 (alpha_old(udd_bp1_0_3047_9665[cond2])),color=color_uddbp1, ls = '-.', lw=lw_uddbp1)
p6 = a.fill(np.NaN, np.NaN, color=color_uddbp1, alpha=0.3)

a.fill_between(np.array(lle_bp2_0_3047_9665[cond3]['phi0'], dtype=float),
                 (alpha_old(lle_bp2_0_3047_9627[cond3])), 
                 (alpha_old(lle_bp2_0_3047_9703[cond3])),
                 color=color_llebp2, alpha=0.3)
p7 = a.plot(np.array(lle_bp2_0_3047_9665[cond3]['phi0'], dtype=float),
                 (alpha_old(lle_bp2_0_3047_9665[cond3])),color=color_llebp2, ls = '-.', lw=lw_uddbp1)
p8 = a.fill(np.NaN, np.NaN, color=color_llebp2, alpha=0.3)

a.fill_between(np.array(udd_bp2_0_3047_9665[cond4]['phi0'], dtype=float),
                 (alpha_old(udd_bp2_0_3047_9627[cond4])), 
                 (alpha_old(udd_bp2_0_3047_9703[cond4])),
                 color=color_uddbp2, alpha=0.3)
p9 = a.plot(np.array(udd_bp2_0_3047_9665[cond2]['phi0'], dtype=float),
                 (alpha_old(udd_bp2_0_3047_9665[cond4])),color=color_uddbp2, ls = '-.', lw=lw_uddbp1)
p10 = a.fill(np.NaN, np.NaN, color=color_uddbp2, alpha=0.3)

a.set_xlabel(r'$\phi_0$ (GeV)')
a.set_ylabel(r'$\alpha\equiv1-\frac{A_6(\phi_0)^2}{20m_\phi(\phi_0)^2}$')
a.set_xscale('log')
#a.set_yscale('log')

a.legend([(p2[0], p1[0]), (p3[0], p4[0]), (p5[0], p6[0]), (p7[0], p8[0]), (p9[0], p10[0]), ], ['tree','lle bp1','udd bp1','lle bp2','udd bp2'], fontsize = 25)

a.axvspan(3e16, 4e16, facecolor = 'white')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
a.set_xlim(1e14,3.5e16)
a.set_ylim(-0.01,0.04)
plt.grid()
define_plot_resolution()

plt.savefig('save_plots/ft_old.pdf', bbox_inches='tight')
print('ft_old.pdf saved')
#plt.show()


f, a = plt.subplots()
a.fill_between(np.array(tree_0_3047_9665['phi0'], dtype=float),
                 [float(1-mp.mpf(x)/4) for x in tree_0_3047_9627['alpha']], 
                 [float(1-mp.mpf(x)/4) for x in tree_0_3047_9703['alpha']],
                 color=color_tree, alpha=0.3)
p1 = a.plot(np.array(tree_0_3047_9665['phi0'], dtype=float),
                 [float(1-mp.mpf(x)) for x in tree_0_3047_9665['alpha']],color = color_tree)
p2 = a.fill(np.NaN, np.NaN, color_tree, alpha=0.3)
a.fill_between(np.array(lle_bp1_0_3047_9665[cond1]['phi0'], dtype=float),
                 [float(1-mp.mpf(x)/4) for x in lle_bp1_0_3047_9627[cond1]['alpha']], 
                 [float(1-mp.mpf(x)/4) for x in lle_bp1_0_3047_9703[cond1]['alpha']],
                 color=color_llebp1, alpha=0.3)
p3 = a.plot(np.array(lle_bp1_0_3047_9665[cond1]['phi0'], dtype=float),
                 [float(1-mp.mpf(x)/4) for x in lle_bp1_0_3047_9665[cond1]['alpha']],color=color_llebp1, ls = ':')
p4 = a.fill(np.NaN, np.NaN, color=color_llebp1, alpha=0.3)

a.fill_between(np.array(udd_bp1_0_3047_9665[cond2]['phi0'], dtype=float),
                 [float(1-mp.mpf(x)/4) for x in udd_bp1_0_3047_9627[cond2]['alpha']], 
                 [float(1-mp.mpf(x)) for x in udd_bp1_0_3047_9703[cond2]['alpha']],
                 color=color_uddbp1, alpha=0.3)
p5 = a.plot(np.array(udd_bp1_0_3047_9665[cond2]['phi0'], dtype=float),
                 [float(1-mp.mpf(x)/4) for x in udd_bp1_0_3047_9665[cond2]['alpha']],color=color_uddbp1, ls = '-.')
p6 = a.fill(np.NaN, np.NaN, color=color_uddbp1, alpha=0.3)

a.fill_between(np.array(lle_bp2_0_3047_9665[cond3]['phi0'], dtype=float),
                 [float(1-mp.mpf(x)/4) for x in lle_bp2_0_3047_9627[cond3]['alpha']], 
                 [float(1-mp.mpf(x)/4) for x in lle_bp2_0_3047_9703[cond3]['alpha']],
                 color=color_llebp2, alpha=0.3)
p7 = a.plot(np.array(lle_bp2_0_3047_9665[cond3]['phi0'], dtype=float),
                 [float(1-mp.mpf(x)/4) for x in lle_bp2_0_3047_9665[cond3]['alpha']],color=color_llebp2, ls = ':')
p8 = a.fill(np.NaN, np.NaN, color=color_llebp2, alpha=0.3)

a.fill_between(np.array(udd_bp2_0_3047_9665[cond4]['phi0'], dtype=float),
                 [float(1-mp.mpf(x)/4) for x in udd_bp2_0_3047_9627[cond4]['alpha']], 
                 [float(1-mp.mpf(x)/4) for x in udd_bp2_0_3047_9703[cond4]['alpha']],
                 color=color_uddbp2, alpha=0.3)
p9 = a.plot(np.array(udd_bp2_0_3047_9665[cond4]['phi0'], dtype=float),
                 [float(1-mp.mpf(x)/4) for x in udd_bp2_0_3047_9665[cond4]['alpha']],color=color_uddbp2, ls = '-.')
p10 = a.fill(np.NaN, np.NaN, color=color_uddbp2, alpha=0.3)

a.set_xlabel(r'$\phi_0$ (GeV)')
a.set_ylabel(r'$\alpha^{\mathrm{(loop)}}$')
a.set_xscale('log')
a.set_yscale('log')

a.legend([(p2[0], p1[0]), (p3[0], p4[0]), (p5[0], p6[0]), (p7[0], p8[0]), (p9[0], p10[0]), ], ['tree','lle bp1','udd bp1','lle bp2','udd bp2'], fontsize = 25)

a.axvspan(3e16, 4e16, facecolor = 'white')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
a.set_xlim(1e14,3.5e16)
a.set_ylim(1e-23,1e-13)
plt.grid()
define_plot_resolution()

plt.savefig('save_plots/ft_new.pdf', bbox_inches='tight')
print('ft_new.pdf saved')
#plt.show()

####

plt.figure()

plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], lle_bp1_0_3047_9665[cond1]['mphi_phi0']-tree_0_3047_9665_high_llebp1[cond1]['mphi_phi0'], label = r'$\Delta_{lle}^{BP_1}$', ls = ls_llebp1, color = color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], udd_bp1_0_3047_9665[cond2]['mphi_phi0']-tree_0_3047_9665_high_uddbp1[cond2]['mphi_phi0'], label = r'$\Delta_{udd}^{BP_1}$', ls = ls_uddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], lle_bp2_0_3047_9665[cond3]['mphi_phi0']-tree_0_3047_9665_high_llebp2[cond3]['mphi_phi0'], label = r'$\Delta_{lle}^{BP_2}$', ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], udd_bp2_0_3047_9665[cond4]['mphi_phi0']-tree_0_3047_9665_high_uddbp2[cond4]['mphi_phi0'], label = r'$\Delta_{udd}^{BP_2}$', ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)
plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], abs(lle_bp1_0_3047_9703[cond1]['mphi_phi0']-lle_bp1_0_3047_9627[cond1]['mphi_phi0'])/2, label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color =color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], abs(udd_bp1_0_3047_9703[cond2]['mphi_phi0']-udd_bp1_0_3047_9627[cond2]['mphi_phi0'])/2, label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], abs(lle_bp2_0_3047_9703[cond3]['mphi_phi0']-lle_bp2_0_3047_9627[cond3]['mphi_phi0'])/2, label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], abs(udd_bp2_0_3047_9703[cond4]['mphi_phi0']-udd_bp2_0_3047_9627[cond4]['mphi_phi0'])/2, label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)


plt.loglog()
plt.legend(fontsize = 20, ncol=2)
plt.axvspan(3e16, 4e16, facecolor = 'white')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.grid()
plt.xlim(1e14,3.5e16)
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'Error on $m_\phi(\phi_0)$ (GeV)')

define_plot_resolution()
plt.savefig('save_plots/error_mphi_phi0_phi0.pdf', bbox_inches='tight')
print('error_mphi_phi0_phi0.pdf saved')


plt.figure()


plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], (lle_bp1_0_3047_9665[cond1]['A6_phi0']-tree_0_3047_9665_high_llebp1[cond1]['A6_phi0']), label = r'$\Delta_{lle}^{BP_1}$', ls = ls_llebp1, color = color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], (udd_bp1_0_3047_9665[cond2]['A6_phi0']-tree_0_3047_9665_high_uddbp1[cond2]['A6_phi0']), label = r'$\Delta_{udd}^{BP_1}$', ls = ls_uddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], (lle_bp2_0_3047_9665[cond3]['A6_phi0']-tree_0_3047_9665_high_llebp2[cond3]['A6_phi0']), label = r'$\Delta_{lle}^{BP_2}$', ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], (udd_bp2_0_3047_9665[cond4]['A6_phi0']-tree_0_3047_9665_high_uddbp2[cond4]['A6_phi0']), label = r'$\Delta_{udd}^{BP_2}$', ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)
plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], abs(lle_bp1_0_3047_9703[cond1]['A6_phi0']-lle_bp1_0_3047_9627[cond1]['A6_phi0'])/2, label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color =color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], abs(udd_bp1_0_3047_9703[cond2]['A6_phi0']-udd_bp1_0_3047_9627[cond2]['A6_phi0'])/2, label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], abs(lle_bp2_0_3047_9703[cond3]['A6_phi0']-lle_bp2_0_3047_9627[cond3]['A6_phi0'])/2, label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], abs(udd_bp2_0_3047_9703[cond4]['A6_phi0']-udd_bp2_0_3047_9627[cond4]['A6_phi0'])/2, label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)

plt.loglog()
plt.legend(fontsize = 20, ncol=2)
plt.axvspan(3e16, 4e16, facecolor = 'white')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.xlim(1e14,3.5e16)
plt.grid()
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'Error on $A_6(\phi_0)$ (GeV)')

define_plot_resolution()
plt.savefig('save_plots/error_A6_phi0_phi0.pdf', bbox_inches='tight')
print('error_A6_phi0_phi0.pdf saved')

plt.figure()
#fig, ax = plt.subplots()

plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], (lle_bp1_0_3047_9665[cond1]['A6_phi0']-tree_0_3047_9665_high_llebp1[cond1]['A6_phi0']), label = r'$\Delta_{lle}^{BP_1}$', ls = ls_llebp1, color = color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], (udd_bp1_0_3047_9665[cond2]['A6_phi0']-tree_0_3047_9665_high_uddbp1[cond2]['A6_phi0']), label = r'$\Delta_{udd}^{BP_1}$', ls = ls_uddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], (lle_bp2_0_3047_9665[cond3]['A6_phi0']-tree_0_3047_9665_high_llebp2[cond3]['A6_phi0']), label = r'$\Delta_{lle}^{BP_2}$',ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], (udd_bp2_0_3047_9665[cond4]['A6_phi0']-tree_0_3047_9665_high_uddbp2[cond4]['A6_phi0']), label = r'$\Delta_{udd}^{BP_2}$', ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)
plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], (lle_bp1_0_3047_9703[cond1]['A6_phi0']-lle_bp1_0_3047_9627[cond1]['A6_phi0'])/2,  label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color = color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], (udd_bp1_0_3047_9703[cond2]['A6_phi0']-udd_bp1_0_3047_9627[cond2]['A6_phi0'])/2, label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], (lle_bp2_0_3047_9703[cond3]['A6_phi0']-lle_bp2_0_3047_9627[cond3]['A6_phi0'])/2, label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], (udd_bp2_0_3047_9703[cond4]['A6_phi0']-udd_bp2_0_3047_9627[cond4]['A6_phi0'])/2, label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)
#plt.yscale('symlog')
#plt.xscale('log')
plt.loglog()
plt.legend(fontsize = 20, ncol=2)
plt.axvspan(3e16, 4e16, facecolor = 'white')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.grid()
plt.xlim(1e14,3.5e16)
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'Error on $A_6(\phi_0)$ (GeV)')

define_plot_resolution()

sub_axes = plt.axes([.612, .22, .32, .35])
sub_axes.plot(lle_bp1_0_3047_9665[cond1]['phi0'], (lle_bp1_0_3047_9665[cond1]['A6_phi0']-tree_0_3047_9665_high_llebp1[cond1]['A6_phi0']), label = r'$\Delta_{lle}^{BP_1}$',ls = ls_llebp1, color = color_llebp1, lw = lw_llebp1)
sub_axes.plot(udd_bp1_0_3047_9665[cond2]['phi0'], (udd_bp1_0_3047_9665[cond2]['A6_phi0']-tree_0_3047_9665_high_uddbp1[cond2]['A6_phi0']), label = r'$\Delta_{udd}^{BP_1}$',ls = ls_uddbp1, color = color_uddbp1, lw = lw_uddbp1)
sub_axes.plot(lle_bp2_0_3047_9665[cond3]['phi0'], (lle_bp2_0_3047_9665[cond3]['A6_phi0']-tree_0_3047_9665_high_llebp2[cond3]['A6_phi0']), label = r'$\Delta_{lle}^{BP_2}$',ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
sub_axes.plot(udd_bp2_0_3047_9665[cond4]['phi0'], (udd_bp2_0_3047_9665[cond4]['A6_phi0']-tree_0_3047_9665_high_uddbp2[cond4]['A6_phi0']), label = r'$\Delta_{udd}^{BP_2}$',ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)
sub_axes.plot(lle_bp1_0_3047_9665[cond1]['phi0'], (lle_bp1_0_3047_9703[cond1]['A6_phi0']-lle_bp1_0_3047_9627[cond1]['A6_phi0'])/2, label = r'$\sigma_{n_s,lle}^{BP_1}$',ls = ls_sllebp1, color = color_llebp1, lw = lw_llebp1)
sub_axes.plot(udd_bp1_0_3047_9665[cond2]['phi0'], (udd_bp1_0_3047_9703[cond2]['A6_phi0']-udd_bp1_0_3047_9627[cond2]['A6_phi0'])/2, label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
sub_axes.plot(lle_bp2_0_3047_9665[cond3]['phi0'], (lle_bp2_0_3047_9703[cond3]['A6_phi0']-lle_bp2_0_3047_9627[cond3]['A6_phi0'])/2, label = r'$\sigma_{n_s,lle}^{BP_2}$',ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2)
sub_axes.plot(udd_bp2_0_3047_9665[cond4]['phi0'], (udd_bp2_0_3047_9703[cond4]['A6_phi0']-udd_bp2_0_3047_9627[cond4]['A6_phi0'])/2, label = r'$\sigma_{n_s,udd}^{BP_2}$',ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)
#plt.xscale('log')
plt.xlim(1e14,.5e15)
plt.ylim(-30,320)
#plt.xlabel(r'$\phi_0$ (GeV)')
#plt.ylabel(r'Error on $A_6(\phi_0)$ (GeV)')
plt.locator_params(axis='y', nbins=4)
plt.xticks([1e14,3e14,5e14], [r'$10^{14}$', r'$3\times10^{14}$', r'$5\times10^{14}$'])
plt.grid()

define_plot_resolution()
#plt.show()
plt.tight_layout()
plt.savefig('save_plots/error_A6_phi0_phi0.pdf', bbox_inches='tight')
print('error_A6_phi0_phi0.pdf saved')

plt.figure()

plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], lle_bp1_0_3047_9665[cond1]['lambda6_phi0']-tree_0_3047_9665_high_llebp1[cond1]['lambda6_phi0'], label = r'$\Delta_{lle}^{BP_1}$', ls = ls_llebp1, color = color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], udd_bp1_0_3047_9665[cond2]['lambda6_phi0']-tree_0_3047_9665_high_uddbp1[cond2]['lambda6_phi0'], label = r'$\Delta_{udd}^{BP_1}$', ls = ls_uddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], lle_bp2_0_3047_9665[cond3]['lambda6_phi0']-tree_0_3047_9665_high_llebp2[cond3]['lambda6_phi0'], label = r'$\Delta_{lle}^{BP_2}$', ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], udd_bp2_0_3047_9665[cond4]['lambda6_phi0']-tree_0_3047_9665_high_uddbp2[cond4]['lambda6_phi0'], label = r'$\Delta_{udd}^{BP_2}$', ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)
plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], abs(lle_bp1_0_3047_9703[cond1]['lambda6_phi0']-lle_bp1_0_3047_9627[cond1]['lambda6_phi0'])/2, label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color =color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], abs(udd_bp1_0_3047_9703[cond2]['lambda6_phi0']-udd_bp1_0_3047_9627[cond2]['lambda6_phi0'])/2, label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], abs(lle_bp2_0_3047_9703[cond3]['lambda6_phi0']-lle_bp2_0_3047_9627[cond3]['lambda6_phi0'])/2, label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], abs(udd_bp2_0_3047_9703[cond4]['lambda6_phi0']-udd_bp2_0_3047_9627[cond4]['lambda6_phi0'])/2, label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)

plt.loglog()
plt.legend(fontsize = 20, ncol=2)
plt.axvspan(3e16, 4e16, facecolor = 'white')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.xlim(1e14,3.5e16)
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'Error on $\lambda_6(\phi_0)$')
plt.grid()
#
define_plot_resolution()
plt.savefig('save_plots/error_lambda6_phi0_phi0.pdf', bbox_inches='tight')
print('error_lambda6_phi0_phi0.pdf saved')


plt.figure()

plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], lle_bp1_0_3047_9665[cond1]['mphi_2']-tree_0_3047_9665_high_llebp1[cond1]['mphi_2'], label = r'$\Delta_{lle}^{BP_1}$', ls = ls_llebp1, color = color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], udd_bp1_0_3047_9665[cond2]['mphi_2']-tree_0_3047_9665_high_uddbp1[cond2]['mphi_2'], label = r'$\Delta_{udd}^{BP_1}$', ls = ls_uddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], lle_bp2_0_3047_9665[cond3]['mphi_2']-tree_0_3047_9665_high_llebp2[cond3]['mphi_2'], label = r'$\Delta_{lle}^{BP_2}$', ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], udd_bp2_0_3047_9665[cond4]['mphi_2']-tree_0_3047_9665_high_uddbp2[cond4]['mphi_2'], label = r'$\Delta_{udd}^{BP_2}$', ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)
plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], abs(lle_bp1_0_3047_9703[cond1]['mphi_2']-lle_bp1_0_3047_9627[cond1]['mphi_2'])/2, label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color =color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], abs(udd_bp1_0_3047_9703[cond2]['mphi_2']-udd_bp1_0_3047_9627[cond2]['mphi_2'])/2, label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], abs(lle_bp2_0_3047_9703[cond3]['mphi_2']-lle_bp2_0_3047_9627[cond3]['mphi_2'])/2, label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], abs(udd_bp2_0_3047_9703[cond4]['mphi_2']-udd_bp2_0_3047_9627[cond4]['mphi_2'])/2, label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)

plt.loglog()
plt.legend(fontsize = 20, ncol=2)
plt.axvspan(3e16, 4e16, facecolor = 'white')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.xlim(1e14,3.5e16)
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'Error on $m_\phi($2 TeV$)$ (GeV)')
plt.grid()

define_plot_resolution()
plt.savefig('save_plots/error_mphi_phi0_2000.pdf', bbox_inches='tight')
print('error_mphi_phi0_2000.pdf saved')


plt.figure()

plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], (lle_bp1_0_3047_9665[cond1]['A6_2']-tree_0_3047_9665_high_llebp1[cond1]['A6_2']), label = r'$\Delta_{lle}^{BP_1}$', ls = ls_llebp1, color = color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], (udd_bp1_0_3047_9665[cond2]['A6_2']-tree_0_3047_9665_high_uddbp1[cond2]['A6_2']), label = r'$\Delta_{udd}^{BP_1}$', ls = ls_uddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], (lle_bp2_0_3047_9665[cond3]['A6_2']-tree_0_3047_9665_high_llebp2[cond3]['A6_2']), label = r'$\Delta_{lle}^{BP_2}$', ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], (udd_bp2_0_3047_9665[cond4]['A6_2']-tree_0_3047_9665_high_uddbp2[cond4]['A6_2']), label = r'$\Delta_{udd}^{BP_2}$', ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)
plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], abs(lle_bp1_0_3047_9703[cond1]['A6_2']-lle_bp1_0_3047_9627[cond1]['A6_2'])/2, label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color =color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], abs(udd_bp1_0_3047_9703[cond2]['A6_2']-udd_bp1_0_3047_9627[cond2]['A6_2'])/2, label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], abs(lle_bp2_0_3047_9703[cond3]['A6_2']-lle_bp2_0_3047_9627[cond3]['A6_2'])/2, label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], abs(udd_bp2_0_3047_9703[cond4]['A6_2']-udd_bp2_0_3047_9627[cond4]['A6_2'])/2, label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)

plt.loglog()
plt.legend(fontsize = 20, ncol=2)
plt.axvspan(3e16, 4e16, facecolor = 'white')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.xlim(1e14,3.5e16)
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'Error on $A_6($2 TeV$)$ (GeV)')
plt.grid()

define_plot_resolution()
plt.savefig('save_plots/error_A6_phi0_2000.pdf', bbox_inches='tight')
print('error_A6_phi0_2000.pdf saved')


plt.figure()

plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], lle_bp1_0_3047_9665[cond1]['lambda6_2']-tree_0_3047_9665_high_llebp1[cond1]['lambda6_2'], label = r'$\Delta_{lle}^{BP_1}$', ls = ls_llebp1, color = color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], udd_bp1_0_3047_9665[cond2]['lambda6_2']-tree_0_3047_9665_high_uddbp1[cond2]['lambda6_2'], label = r'$\Delta_{udd}^{BP_1}$', ls = ls_uddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], lle_bp2_0_3047_9665[cond3]['lambda6_2']-tree_0_3047_9665_high_llebp2[cond3]['lambda6_2'], label = r'$\Delta_{lle}^{BP_2}$', ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], udd_bp2_0_3047_9665[cond4]['lambda6_2']-tree_0_3047_9665_high_uddbp2[cond4]['lambda6_2'], label = r'$\Delta_{udd}^{BP_2}$', ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)
plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], abs(lle_bp1_0_3047_9703[cond1]['lambda6_2']-lle_bp1_0_3047_9627[cond1]['lambda6_2'])/2, label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color =color_llebp1, lw = lw_llebp1)
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], abs(udd_bp1_0_3047_9703[cond2]['lambda6_2']-udd_bp1_0_3047_9627[cond2]['lambda6_2'])/2, label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], abs(lle_bp2_0_3047_9703[cond3]['lambda6_2']-lle_bp2_0_3047_9627[cond3]['lambda6_2'])/2, label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], abs(udd_bp2_0_3047_9703[cond4]['lambda6_2']-udd_bp2_0_3047_9627[cond4]['lambda6_2'])/2, label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)

plt.loglog()
plt.legend(fontsize = 20, ncol=2)
plt.axvspan(3e16, 4e16, facecolor = 'white')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.xlim(1e14,3.5e16)
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'Error on $\lambda_6($2 TeV$)$')
plt.grid()

define_plot_resolution()
plt.savefig('save_plots/error_lambda6_phi0_2000.pdf', bbox_inches='tight')
print('error_lambda6_phi0_2000.pdf saved')
#######

plt.figure()

plt.plot(lle_bp1_0_3047_9665[cond1]['phi0'], (lle_bp1_0_3047_9665[cond1]['mphi_phi0']-tree_0_3047_9665_high_llebp1[cond1]['mphi_phi0'])/tree_0_3047_9665_high_llebp1[cond1]['mphi_phi0'], label = r'$\Delta_{lle}^{BP_1}$', ls = ':', color = 'blue')
plt.plot(udd_bp1_0_3047_9665[cond2]['phi0'], (udd_bp1_0_3047_9665[cond2]['mphi_phi0']-tree_0_3047_9665_high_uddbp1[cond2]['mphi_phi0'])/tree_0_3047_9665_high_uddbp1[cond2]['mphi_phi0'], label = r'$\Delta_{udd}^{BP_1}$', ls = ':', color = 'red')
plt.plot(lle_bp2_0_3047_9665[cond3]['phi0'], (lle_bp2_0_3047_9665[cond3]['mphi_phi0']-tree_0_3047_9665_high_llebp2[cond3]['mphi_phi0'])/tree_0_3047_9665_high_llebp2[cond3]['mphi_phi0'], label = r'$\Delta_{lle}^{BP_2}$', ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(udd_bp2_0_3047_9665[cond4]['phi0'], (udd_bp2_0_3047_9665[cond4]['mphi_phi0']-tree_0_3047_9665_high_uddbp2[cond4]['mphi_phi0'])/tree_0_3047_9665_high_uddbp2[cond4]['mphi_phi0'], label = r'$\Delta_{udd}^{BP_2}$', ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)
plt.plot(pd_to_array(tree_0_3047_9665, 'phi0'), abs(pd_to_array(tree_0_3047_9703, 'mphi_phi0')-pd_to_array(tree_0_3047_9627, 'mphi_phi0'))/(2*pd_to_array(tree_0_3047_9665, 'mphi_phi0')), label = r'$\sigma_{n_s}$', ls = '--', color = 'black')
plt.plot(pd_to_array(tree_0_3047_9665, 'phi0'), abs(pd_to_array(tree_0_3061_9665, 'mphi_phi0')-pd_to_array(tree_0_3033_9665, 'mphi_phi0'))/(2*pd_to_array(tree_0_3047_9665, 'mphi_phi0')), label = r'$\sigma_{A_s}$', ls = '-', color = 'black')
plt.plot(pd_to_array(tree_0_3047_9665, 'phi0'), abs(pd_to_array(tree_0_3047_9665, 'mphi_phi0')-pd_to_array(tree_10_3047_9665, 'mphi_phi0'))/(10*pd_to_array(tree_0_3047_9665, 'mphi_phi0')), label = r'$\sigma_{lnRrad\pm1}$', ls = '-.', color = 'black')

plt.loglog()
plt.legend(fontsize = 20, ncol=2)
plt.axvspan(3e16, 4e16, facecolor = 'white')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.xlim(1e14,3.5e16)
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'Relative error on $m_\phi(\phi_0)$ (GeV)')
plt.grid()

define_plot_resolution()
plt.savefig('save_plots/complete_error_mphi_phi0_phi0.pdf', bbox_inches='tight')
print('complete_error_mphi_phi0_phi0.pdf saved')


plt.figure() # to zoom on a particular zone

#plt.plot(tree_0_3047_9665['A6_phi0'], tree_0_3047_9665['mphi_phi0'], ls = '-', color = color_tree)
#plt.plot(tree_0_3047_9627['A6_phi0'], tree_0_3047_9627['mphi_phi0'], label = r'tree', ls = '-', color = color_tree)
#plt.plot(tree_0_3047_9703['A6_phi0'], tree_0_3047_9703['mphi_phi0'], ls = '-', color = color_tree)

plt.plot(tree_0_3047_9665_high_uddbp2['A6_phi0'], tree_0_3047_9665_high_uddbp2['mphi_phi0'], ls = '-', color = color_tree)
plt.plot(tree_0_3047_9665_high_uddbp2['A6_phi0'], tree_0_3047_9665_high_uddbp2['mphi_phi0'], label = r'tree', ls = '-', color = color_tree)
plt.plot(tree_0_3047_9665_high_uddbp2['A6_phi0'], tree_0_3047_9665_high_uddbp2['mphi_phi0'], ls = '-', color = color_tree)

plt.plot(lle_bp1_0_3047_9665[cond1]['A6_phi0'], lle_bp1_0_3047_9665[cond1]['mphi_phi0'], ls = ':', color = color_llebp1)
plt.plot(lle_bp1_0_3047_9627[cond1]['A6_phi0'], lle_bp1_0_3047_9627[cond1]['mphi_phi0'], label = r'lle BP$_1$', ls = '-', color = color_llebp1)
plt.plot(lle_bp1_0_3047_9703[cond1]['A6_phi0'], lle_bp1_0_3047_9703[cond1]['mphi_phi0'], ls = '-.', color = color_llebp1)

plt.plot(udd_bp1_0_3047_9665[cond2]['A6_phi0'], udd_bp1_0_3047_9665[cond2]['mphi_phi0'], ls = ':', color = color_uddbp1)
plt.plot(udd_bp1_0_3047_9627[cond2]['A6_phi0'], udd_bp1_0_3047_9627[cond2]['mphi_phi0'], label = r'udd BP$_1$', ls = '-', color = color_uddbp1)
plt.plot(udd_bp1_0_3047_9703[cond2]['A6_phi0'], udd_bp1_0_3047_9703[cond2]['mphi_phi0'], ls = '-.', color = color_uddbp1)

plt.plot(lle_bp2_0_3047_9665[cond3]['A6_phi0'], lle_bp2_0_3047_9665[cond3]['mphi_phi0'], ls = ':', color = color_llebp2)
plt.plot(lle_bp2_0_3047_9627[cond3]['A6_phi0'], lle_bp2_0_3047_9627[cond3]['mphi_phi0'], label = r'lle BP$_2$', ls = '-', color = color_llebp2)
plt.plot(lle_bp2_0_3047_9703[cond3]['A6_phi0'], lle_bp2_0_3047_9703[cond3]['mphi_phi0'], ls = '-.', color = color_llebp2)

plt.plot(udd_bp2_0_3047_9665[cond4]['A6_phi0'], udd_bp2_0_3047_9665[cond4]['mphi_phi0'], ls = ':', color = color_uddbp2)
plt.plot(udd_bp2_0_3047_9627[cond4]['A6_phi0'], udd_bp2_0_3047_9627[cond4]['mphi_phi0'], label = r'udd BP$_2$', ls = '-', color = color_uddbp2)
plt.plot(udd_bp2_0_3047_9703[cond4]['A6_phi0'], udd_bp2_0_3047_9703[cond4]['mphi_phi0'], ls = '-.', color = color_uddbp2)

plt.loglog()
plt.legend(fontsize = 20)
plt.xlabel(r'$A_6$ (GeV)')
plt.ylabel(r'$m_\phi$ (GeV)')
plt.grid()

define_plot_resolution()
plt.savefig('save_plots/A6_mphi.pdf', bbox_inches='tight')
print('A6_mphi.pdf saved')

#plt.figure()
#plt.plot(pd_to_array(tree_0_3047_9665, 'phi0'), pd_to_array(tree_0_3047_9665, 'A6_phi0')/1000,'red')
#plt.fill_between(pd_to_array(tree_0_3047_9665, 'phi0'), pd_to_array(tree_0_3047_9627, 'A6_phi0')/1000, pd_to_array(tree_0_3047_9703, 'A6_phi0')/1000, color='red', alpha=0.2, label=r'$lnR_{rad} = 0$ ; $\overline{n_s}\pm \sigma_{n_s}$')
#plt.plot(pd_to_array(lle_bp1_0_3047_9665, 'phi0'), pd_to_array(lle_bp1_0_3047_9665, 'A6_phi0')/1000,'yellow')
#plt.fill_between(pd_to_array(lle_bp1_0_3047_9665, 'phi0'), pd_to_array(lle_bp1_0_3047_9627, 'A6_phi0')/1000, pd_to_array(lle_bp1_0_3047_9703, 'A6_phi0')/1000, color='yellow', alpha=0.2, label=r'lle bp1')
#plt.plot(pd_to_array(udd_bp1_0_3047_9665, 'phi0'), pd_to_array(udd_bp1_0_3047_9665, 'A6_phi0')/1000,'green')
#plt.fill_between(pd_to_array(udd_bp1_0_3047_9665, 'phi0'), pd_to_array(udd_bp1_0_3047_9627, 'A6_phi0')/1000, pd_to_array(udd_bp1_0_3047_9703, 'A6_phi0')/1000, color='green', alpha=0.2, label=r'udd bp1')
#plt.plot(pd_to_array(lle_bp2_0_3047_9665, 'phi0'), pd_to_array(lle_bp2_0_3047_9665, 'A6_phi0')/1000,'blue')
#plt.fill_between(pd_to_array(lle_bp2_0_3047_9627, 'phi0'), pd_to_array(lle_bp2_0_3047_9627, 'A6_phi0')/1000, pd_to_array(lle_bp2_0_3047_9703, 'A6_phi0')/1000, color='blue', alpha=0.2, label=r'lle bp2')
#plt.plot(pd_to_array(udd_bp2_0_3047_9665, 'phi0'), pd_to_array(udd_bp2_0_3047_9665, 'A6_phi0')/1000,'pink')
#plt.fill_between(pd_to_array(udd_bp2_0_3047_9627, 'phi0'), pd_to_array(udd_bp2_0_3047_9627, 'A6_phi0')/1000, pd_to_array(udd_bp2_0_3047_9703, 'A6_phi0')/1000, color='pink', alpha=0.2, label=r'udd bp2')
#plt.xlabel(r'$\phi_0$ (GeV)')
#plt.ylabel(r'$A_{6}$ (TeV)')
#plt.xlim(min(pd_to_array(tree_10_3047_9627, 'phi0')), 8e15)
#plt.legend(fontsize=25)
#plt.ylim(1.8e-1,8e2)
#plt.grid()
#define_plot_resolution()
#
#plt.savefig('save_plots/A6_compare_lin.pdf', bbox_inches='tight')
#print('A6_compare_lin.pdf saved')

from scipy.interpolate import interp1d
def interp(inflaton, bp, phi0or2 = 'phi0'):

    if inflaton == 'lle' and bp == 1: file_rge, file_rge_0953, file_tree = lle_bp1_0_3047_9665, lle_bp1_0_3047_9627, tree_0_3047_9665_high_llebp1
    elif inflaton == 'udd' and bp == 1: file_rge, file_rge_0953, file_tree = udd_bp1_0_3047_9665, udd_bp1_0_3047_9627, tree_0_3047_9665_high_uddbp1
    elif inflaton == 'lle' and bp == 2: file_rge, file_rge_0953, file_tree = lle_bp2_0_3047_9665, lle_bp2_0_3047_9627, tree_0_3047_9665_high_llebp2
    elif inflaton == 'udd' and bp == 2: file_rge, file_rge_0953, file_tree = udd_bp2_0_3047_9665, udd_bp2_0_3047_9627, tree_0_3047_9665_high_uddbp2

    cond = file_rge['yep_or_nop']=='yep!'
    min_interp = max([min([float(x) for x in file_rge[cond]['A6_'+phi0or2]]),
                      min([float(x) for x in file_rge_0953[cond]['A6_'+phi0or2]]),
                      min([float(x) for x in file_tree[cond]['A6_'+phi0or2]])])
                      
    max_interp = min([max([float(x) for x in file_rge[cond]['A6_'+phi0or2]]),
                      max([float(x) for x in file_rge_0953[cond]['A6_'+phi0or2]]),
                      max([float(x) for x in file_tree[cond]['A6_'+phi0or2]])])
                      
    onto_interp = np.logspace(np.log10(min_interp), np.log10(max_interp), 5000,)[1:-1]
    #onto_interp = np.logspace(np.log10(float(file_rge[cond]['A6_'+phi0or2].iloc[0])), np.log10(float(file_rge[cond]['A6_'+phi0or2].iloc[-1])),5000,)[30:-399]
    #print(min(onto_interp), max(onto_interp))
    #print(min(np.array([float(x) for x in file_rge[cond]['A6_'+phi0or2]])), max(np.array([float(x) for x in file_rge[cond]['A6_'+phi0or2]])))
    rge_interpd = interp1d(np.array([float(x) for x in file_rge[cond]['A6_'+phi0or2]]), np.array([float(x) for x in file_rge[cond]['mphi_'+phi0or2]]), kind='cubic')
    #print(min(np.array([float(x) for x in file_rge_0953[cond]['A6_'+phi0or2]])), max(np.array([float(x) for x in file_rge_0953[cond]['A6_'+phi0or2]])))
    rge_interpd_0953 = interp1d(np.array([float(x) for x in file_rge_0953[cond]['A6_'+phi0or2]]), np.array([float(x) for x in file_rge_0953[cond]['mphi_'+phi0or2]]), kind='cubic')
    #print(min(np.array([float(x) for x in file_tree[cond]['A6_'+phi0or2]])), max(np.array([float(x) for x in file_tree[cond]['A6_'+phi0or2]])))
    tree_interpd = interp1d(np.array([float(x) for x in file_tree[cond]['A6_'+phi0or2]]), np.array([float(x) for x in file_tree[cond]['mphi_'+phi0or2]]), kind='cubic')
    #print('')
    
    rge_interpd = rge_interpd(onto_interp)
    rge_interpd_0953 = rge_interpd_0953(onto_interp)
    tree_interpd = tree_interpd(onto_interp)
    
    return onto_interp, rge_interpd, tree_interpd, rge_interpd_0953

onto_lle_bp1_0_3047_9665_interp, lle_bp1_0_3047_9665_interp, tree_0_3047_9665_high_llebp1_interp, lle_bp1_0_3047_9627_interp = interp('lle', 1, 'phi0')
onto_udd_bp1_0_3047_9665_interp, udd_bp1_0_3047_9665_interp, tree_0_3047_9665_high_uddbp1_interp, udd_bp1_0_3047_9627_interp = interp('udd', 1, 'phi0')
onto_lle_bp2_0_3047_9665_interp, lle_bp2_0_3047_9665_interp, tree_0_3047_9665_high_llebp2_interp, lle_bp2_0_3047_9627_interp = interp('lle', 2, 'phi0')
onto_udd_bp2_0_3047_9665_interp, udd_bp2_0_3047_9665_interp, tree_0_3047_9665_high_uddbp2_interp, udd_bp2_0_3047_9627_interp = interp('udd', 2, 'phi0')


plt.figure()

plt.plot(onto_lle_bp1_0_3047_9665_interp, abs(lle_bp1_0_3047_9665_interp - tree_0_3047_9665_high_llebp1_interp), label = r'$\Delta_{lle}^{BP_1}$', ls = ls_llebp1, color = color_llebp1, lw = lw_llebp1)
plt.plot(onto_udd_bp1_0_3047_9665_interp, abs(udd_bp1_0_3047_9665_interp - tree_0_3047_9665_high_uddbp1_interp), label = r'$\Delta_{udd}^{BP_1}$', ls = ls_uddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(onto_lle_bp2_0_3047_9665_interp, abs(lle_bp2_0_3047_9665_interp - tree_0_3047_9665_high_llebp2_interp), label = r'$\Delta_{lle}^{BP_2}$', ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(onto_udd_bp2_0_3047_9665_interp, abs(udd_bp2_0_3047_9665_interp - tree_0_3047_9665_high_uddbp2_interp), label = r'$\Delta_{udd}^{BP_2}$', ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)

#plt.plot(onto_lle_bp1_0_3047_9665_interp, abs(lle_bp1_0_3047_9665_interp-lle_bp1_0_3047_9627_interp), label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color =color_llebp1, lw = lw_llebp1)
#plt.plot(onto_udd_bp1_0_3047_9665_interp, abs(udd_bp1_0_3047_9665_interp-udd_bp1_0_3047_9627_interp), label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
#plt.plot(onto_lle_bp2_0_3047_9665_interp, abs(lle_bp2_0_3047_9665_interp-lle_bp2_0_3047_9627_interp), label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2
#plt.plot(onto_udd_bp2_0_3047_9665_interp, abs(udd_bp2_0_3047_9665_interp-udd_bp2_0_3047_9627_interp), label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)

plt.plot(onto_lle_bp1_0_3047_9665_interp, 0.000001*lle_bp1_0_3047_9665_interp, label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color =color_llebp1, lw = lw_llebp1)
plt.plot(onto_udd_bp1_0_3047_9665_interp, 0.000001*udd_bp1_0_3047_9665_interp, label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(onto_lle_bp2_0_3047_9665_interp, 0.000001*lle_bp2_0_3047_9665_interp, label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(onto_udd_bp2_0_3047_9665_interp, 0.000001*udd_bp2_0_3047_9665_interp, label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)

plt.loglog()
plt.xlim(1e3,3e7)
plt.grid()
plt.legend(fontsize = 20, ncol=2)
plt.xlabel(r'$A_6$ (GeV)')
plt.ylabel(r'Error on $m_\phi(\phi_0)$ (GeV)')

define_plot_resolution()
plt.savefig('save_plots/error_mphi_A6_phi0.pdf', bbox_inches='tight')
print('error_mphi_A6_phi0.pdf saved')

onto_lle_bp1_0_3047_9665_interp, lle_bp1_0_3047_9665_interp, tree_0_3047_9665_high_llebp1_interp, lle_bp1_0_3047_9627_interp = interp('lle', 1, '2')
onto_udd_bp1_0_3047_9665_interp, udd_bp1_0_3047_9665_interp, tree_0_3047_9665_high_uddbp1_interp, udd_bp1_0_3047_9627_interp = interp('udd', 1, '2')
onto_lle_bp2_0_3047_9665_interp, lle_bp2_0_3047_9665_interp, tree_0_3047_9665_high_llebp2_interp, lle_bp2_0_3047_9627_interp = interp('lle', 2, '2')
onto_udd_bp2_0_3047_9665_interp, udd_bp2_0_3047_9665_interp, tree_0_3047_9665_high_uddbp2_interp, udd_bp2_0_3047_9627_interp = interp('udd', 2, '2')


plt.figure()

plt.plot(onto_lle_bp1_0_3047_9665_interp, abs(lle_bp1_0_3047_9665_interp - tree_0_3047_9665_high_llebp1_interp), label = r'$\Delta_{lle}^{BP_1}$', ls = ls_llebp1, color = color_llebp1, lw = lw_llebp1)
plt.plot(onto_udd_bp1_0_3047_9665_interp, abs(udd_bp1_0_3047_9665_interp - tree_0_3047_9665_high_uddbp1_interp), label = r'$\Delta_{udd}^{BP_1}$', ls = ls_uddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(onto_lle_bp2_0_3047_9665_interp, abs(lle_bp2_0_3047_9665_interp - tree_0_3047_9665_high_llebp2_interp), label = r'$\Delta_{lle}^{BP_2}$', ls = ls_llebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(onto_udd_bp2_0_3047_9665_interp, abs(udd_bp2_0_3047_9665_interp - tree_0_3047_9665_high_uddbp2_interp), label = r'$\Delta_{udd}^{BP_2}$', ls = ls_uddbp2, color = color_uddbp2, lw = lw_uddbp2)

#plt.plot(onto_lle_bp1_0_3047_9665_interp, abs(lle_bp1_0_3047_9665_interp-lle_bp1_0_3047_9627_interp), label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color =color_llebp1, lw = lw_llebp1)
#plt.plot(onto_udd_bp1_0_3047_9665_interp, abs(udd_bp1_0_3047_9665_interp-udd_bp1_0_3047_9627_interp), label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
#plt.plot(onto_lle_bp2_0_3047_9665_interp, abs(lle_bp2_0_3047_9665_interp-lle_bp2_0_3047_9627_interp), label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2
#plt.plot(onto_udd_bp2_0_3047_9665_interp, abs(udd_bp2_0_3047_9665_interp-udd_bp2_0_3047_9627_interp), label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)

plt.plot(onto_lle_bp1_0_3047_9665_interp, 0.000001*lle_bp1_0_3047_9665_interp, label = r'$\sigma_{n_s,lle}^{BP_1}$', ls = ls_sllebp1, color =color_llebp1, lw = lw_llebp1)
plt.plot(onto_udd_bp1_0_3047_9665_interp, 0.000001*udd_bp1_0_3047_9665_interp, label = r'$\sigma_{n_s,udd}^{BP_1}$', ls = ls_suddbp1, color = color_uddbp1, lw = lw_uddbp1)
plt.plot(onto_lle_bp2_0_3047_9665_interp, 0.000001*lle_bp2_0_3047_9665_interp, label = r'$\sigma_{n_s,lle}^{BP_2}$', ls = ls_sllebp2, color = color_llebp2, lw = lw_llebp2)
plt.plot(onto_udd_bp2_0_3047_9665_interp, 0.000001*udd_bp2_0_3047_9665_interp, label = r'$\sigma_{n_s,udd}^{BP_2}$', ls = ls_suddbp2, color = color_uddbp2, lw = lw_uddbp2)

plt.loglog()
plt.xlim(1e2,3e7)
plt.grid()
plt.legend(fontsize = 20, ncol=2)
plt.xlabel(r'$A_6$ (GeV)')
plt.ylabel(r'Error on $m_\phi($2 TeV$)$ (GeV)')

define_plot_resolution()
plt.savefig('save_plots/error_mphi_A6_2000.pdf', bbox_inches='tight')
print('error_mphi_A6_2000.pdf saved')

with extradps(32):
    phigut = mp.mpf('3e16')
    lnMpinGev = mp.mpf('42.334')
    Mp = mp.exp(lnMpinGev)
    pre = mp.mpf('1')/(mp.mpf('8')*mp.pi**2)
    b1, b2, b3 = mp.mpf('33')/mp.mpf('5'), mp.mpf('1'), mp.mpf('-3')
    bp = ''
    fac_corr = mp.sqrt(mp.mpf('2'))
    fac_highscale = mp.mpf('1')

    g1gut = mp.sqrt(mp.mpf('5')/mp.mpf('3'))*mp.mpf('5.45185741e-01')
    g2gut = mp.mpf('6.90473022e-01')
    g3gut = mp.mpf('6.84972506e-01')
    m1gut = mp.mpf('1.36108022e+02')
    m2gut = mp.mpf('1.14286222e+03')
    m3gut = mp.mpf('8.98639714e+02')

    #g1gut = mp.mpf('0.704555660557172')
    #g2gut = mp.mpf('0.690364970285155')
    #g3gut = mp.mpf('0.684720653567032')
    #m1gut = mp.mpf('897.774785812765')
    #m2gut = mp.mpf('1789.66021839594')
    #m3gut = mp.mpf('882.522487633969')



mphi_tree = lambda phi, mphigut : mphigut
lambda6_tree = lambda phi, lambda6gut : lambda6gut
A6_tree = lambda phi, A6gut : A6gut

g1 = lambda phi : g1gut/(mp.sqrt(1-pre*b1*g1gut**2*mp.log(phi/phigut)))
g2 = lambda phi : g2gut/(mp.sqrt(1-pre*b2*g2gut**2*mp.log(phi/phigut)))
g3 = lambda phi : g3gut/(mp.sqrt(1-pre*b3*g3gut**2*mp.log(phi/phigut)))

M1 = lambda phi : m1gut*(g1(phi)/g1gut)**mp.mpf('2')
M2 = lambda phi : m2gut*(g2(phi)/g2gut)**mp.mpf('2')
M3 = lambda phi : m3gut*(g3(phi)/g3gut)**mp.mpf('2')

mphi_lle = lambda phi, mphigut : mp.sqrt(mphigut**2+(m2gut**2-M2(phi)**2)+mp.mpf('1')/11*(m1gut**2-M1(phi)**2))
A6_lle = lambda phi, A6gut : A6gut-mp.mpf('6')*(m2gut-M2(phi))-mp.mpf('6')/11*(m1gut-M1(phi))
lambda6_lle = lambda phi, lambda6gut : lambda6gut*(g2gut/g2(phi))**mp.mpf('6')*(g1gut/g1(phi))**(mp.mpf('6')/11)

mphi_udd = lambda phi, mphigut : mp.sqrt(mphigut**2-mp.mpf('8')/9*(m3gut**2-M3(phi)**2)+mp.mpf('4')/99*(m1gut**2-M1(phi)**2))
A6_udd = lambda phi, A6gut : A6gut+mp.mpf('16')/3*(m3gut-M3(phi))-mp.mpf('8')/33*(m1gut-M1(phi))
lambda6_udd = lambda phi, lambda6gut : lambda6gut*(g3gut/g3(phi))**(mp.mpf('-16')/3)*(g1gut/g1(phi))**(mp.mpf('8')/33)

def V_MSSM(phi, infl_type, mphigut, A6gut, lambda6gut):
    if infl_type == 0 or infl_type == 'tree':
        mphi_func, A6_func, lambda6_func = mphi_tree, A6_tree, lambda6_tree
    elif infl_type == 1 or infl_type == 'lle':
        mphi_func, A6_func, lambda6_func = mphi_lle, A6_lle, lambda6_lle
    elif infl_type == 2 or infl_type == 'udd':
        mphi_func, A6_func, lambda6_func = mphi_udd, A6_udd, lambda6_udd
    else:
        return 'Error: unknown type of inflation'
    lambda6 = lambda6_func(phi, lambda6gut)
    mphi = mphi_func(phi,mphigut)
    A6 = A6_func(phi,A6gut)
    V = mp.mpf('0.5')*mphi**mp.mpf('2')*phi**mp.mpf('2')-fac_corr*fac_highscale*lambda6*A6/(mp.mpf('6')*Mp**mp.mpf('3'))*phi**mp.mpf('6')+fac_highscale**2*lambda6**mp.mpf('2')*phi**mp.mpf('10')/Mp**mp.mpf('6')
    return V

def Vtree(phi, mphi_GUT, A6_GUT, lambda6_GUT):
    return mp.fadd(mp.fsub(mp.fmul(mp.fmul(0.5,mp.power(mphi_GUT,2)),mp.power(phi,2)),mp.fmul(mp.fdiv(mp.sqrt('2')*mp.fmul(lambda6_GUT,A6_GUT),(mp.fmul(6,mp.power(Mp,3)))),mp.power(phi,6))),mp.fmul(mp.fdiv(mp.power(lambda6_GUT,2),mp.power(Mp,6)),mp.power(phi,10)))
def DV(V, phi):
     return (V(phi+1)-V(phi-1))/mp.mpf('2')
def eps1_(V, phi):
    return mp.fdiv(mp.power((mp.fmul((mp.mpf('1')*Mp),mp.fdiv(DV(V, phi),V(phi)))),2),2)
    

def DeltaNstar(file, lnRrad, infl_type):
    if infl_type == 'tree':phistar_list, mphi_list, A6_list, lambda6_list, phi0_list = pd_to_array(file, 'phi*', dtype='mp'), pd_to_array(file, 'mphi_phi0', dtype='mp'), pd_to_array(file, 'A6_phi0', dtype='mp'), pd_to_array(file, 'lambda6_phi0', dtype='mp'), pd_to_array(file, 'phi0', dtype='mp')
    else:
        phistar_list, mphi_list, A6_list, lambda6_list, phi0_list = pd_to_array(file, 'phi*', dtype='mp'), pd_to_array(file, 'mphi_gut', dtype='mp'), pd_to_array(file, 'A6_gut', dtype='mp'), pd_to_array(file, 'lambda6_gut', dtype='mp'), pd_to_array(file, 'phi0', dtype='mp')
    kstar, lnMpcToKappa, HubbleSquareRootOf3OmegaRad, RelatDofRatio = 0.05, 130.282, 7.5437e-63, 1
    Nstar_list = []
    for i, phistar in enumerate(phistar_list):
        phistar = phistar_list[i]
        #print(mphi_list[i], A6_list[i], lambda6_list[i], phi0_list[i])

        V = lambda phi : V_MSSM(phi, infl_type, mphi_list[i], A6_list[i], lambda6_list[i])
        phi_end = mp.findroot(lambda phi: eps1_(V, phi)-1, 0.98*phi0_list[i], tol = 1e-30)

        N0 = log(kstar) - lnMpcToKappa - 0.5*log(HubbleSquareRootOf3OmegaRad) -0.25*log(RelatDofRatio)
        Nstar = lnRrad - N0 - 0.25*mp.log(9/eps1_(V, phistar)/(3-eps1_(V, phi_end))*V(phi_end)/V(phistar))+0.25*mp.log(8*mp.pi**2*mp.mpf('2.0989031673191437e-9'))
        Nstar_list.append(Nstar)
    return Nstar_list

plt.figure()
plt.plot(pd_to_array(lle_bp1_0_3047_9665[cond1], 'phi0'), DeltaNstar(lle_bp1_0_3047_9665[cond1], 0, 'lle'),'orange',ls='--', label=r'Eq.(14) (lle bp1)')
plt.plot(pd_to_array(udd_bp1_0_3047_9665[cond2], 'phi0'), DeltaNstar(udd_bp1_0_3047_9665[cond2], 0, 'udd'),'pink',ls=':', label=r'Eq.(14) (udd bp1)')

V1, V2 = [], []
for i in range(len(pd_to_array(file111, 'phi0'))):
    V1.append(float(Vtree(pd_to_array(file101, 'phi0')[i],pd_to_array(file101, 'mphi_phi0')[i],pd_to_array(file101, 'A6_phi0')[i],pd_to_array(file101, 'lambda6_phi0')[i],)))
    V2.append(float(Vtree(pd_to_array(file111, 'phi0')[i],pd_to_array(file111, 'mphi_phi0')[i],pd_to_array(file111, 'A6_phi0')[i],pd_to_array(file111, 'lambda6_phi0')[i],)))
kstar, lnMpcToKappa, HubbleSquareRootOf3OmegaRad, RelatDofRatio = 0.05, 130.282, 7.5437e-63, 1

N0_SM = log(kstar) - lnMpcToKappa - 0.5*log(HubbleSquareRootOf3OmegaRad) -0.25*log(RelatDofRatio)
RelatDofRatio = 915/4/(427/4)
N0_MSSM = log(kstar) - lnMpcToKappa - 0.5*log(HubbleSquareRootOf3OmegaRad) -0.25*log(RelatDofRatio)

plt.plot(pd_to_array(file111, 'phi0'), -N0_SM-1/4*np.log(9)+1/4*np.log(np.array(V2)/float(Mp)**4),'blue', label=r'$N_0^{SM}-\frac{1}{4}\log{9}+\frac{1}{4}\log{\frac{V(\phi_0)}{M_p^4}}$')
plt.plot(pd_to_array(file111, 'phi0'), -N0_MSSM-1/4*np.log(9)+1/4*np.log(np.array(V2)/float(Mp)**4),'green', label=r'$N_0^{MSSM}-\frac{1}{4}\log{9}+\frac{1}{4}\log{\frac{V(\phi_0)}{M_p^4}}$')

#plt.fill_between(pd_to_array(file101, 'phi0'), [float(x) for x in DeltaNstar(file100, -10)], [float(x) for x in DeltaNstar(file102, -10)], color = 'blue', alpha=0.2, label=r'$lnR_{rad} = -10$ ; $\overline{n_s}\pm \sigma_{n_s}$')
#plt.fill_between(pd_to_array(file111, 'phi0'), [float(x) for x in DeltaNstar(file110, 0)], [float(x) for x in DeltaNstar(file112, 0)], color='red', alpha=0.2, label=r'$lnR_{rad} = 0$ ; $\overline{n_s}\pm \sigma_{n_s}$')
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'$\Delta N*$')
plt.axvline(3e16, color=color_gut, ls =ls_gut, lw=lw_gut)
plt.xlim(1e14,3.5e16)
plt.semilogx()
plt.grid()
define_plot_resolution()
plt.legend(fontsize=25)

plt.savefig('save_plots/Nstar_phi0_tree_reheat.pdf', bbox_inches='tight')
print('Nstar_phi0_tree_reheat.pdf saved')

####

i=33
infl_type = 'lle'
file_i = lle_bp1_0_3047_9665.iloc[i]
file_tree_i = tree_0_3047_9665_high_llebp1.iloc[i]
phi_list = np.linspace(1e12,float(file_tree_i['phi0'])*1.3, 1000)
V1, V2, V3 = [], [], []
for phi in phi_list:
    V1.append(float(V_MSSM(phi, infl_type, file_i['mphi_gut'], file_i['A6_gut'], file_i['lambda6_gut'])))
    V2.append(float(V_MSSM(phi, 'tree', file_tree_i['mphi_phi0'], file_tree_i['A6_phi0'], file_tree_i['lambda6_phi0'])))
    V3.append(float(V_MSSM(phi, infl_type, file_tree_i['mphi_gut'], file_tree_i['A6_gut'], file_tree_i['lambda6_gut'])))

plt.figure()
plt.plot(phi_list, V2, label = r'$V_{\mathrm{tree}}(p_{\mathrm{tree}})$', ls = '-', color = color_llebp1, lw = 2)
plt.plot(phi_list, V1, label = r'$V_{\mathrm{llebp1}}(p_{\mathrm{llebp1}})$', ls = ':', color = color_llebp2, lw = 2)
plt.plot(phi_list, V3, label = r'$V_{\mathrm{llebp1}}(p_{\mathrm{tree}})$', ls = '--', color = color_uddbp1, lw = 2)
plt.axvline(file_tree_i['phi0'], ls=':', lw=4, color='gray')

plt.legend(fontsize = 25)
plt.locator_params(axis='y', nbins=6)
plt.locator_params(axis='x', nbins=6)
plt.text(float(file_tree_i['phi0'])-0.4e14, -0.09e38, r'$\phi_0$', color = 'gray', fontsize=29)
plt.xlabel(r'$\phi$ (GeV)')
plt.ylabel(r'$V(\phi)$ (GeV$^4$)')
plt.xlim(0,float(file_tree_i['phi0'])*1.3)
plt.ylim(bottom = 0)
plt.grid()
define_plot_resolution()

plt.savefig('save_plots/potential.pdf', bbox_inches='tight')
print('potential.pdf saved')


eps1_1, eps1_2, eps1_3, eps1_4 = [], [], [], []
for phi in phi_list:
    eps1_1.append(float(eps1_(lambda phi : V_MSSM(phi, infl_type, file_i['mphi_gut'], file_i['A6_gut'], file_i['lambda6_gut']), phi)))
    eps1_2.append(float(eps1_(lambda phi : V_MSSM(phi, 'tree', file_tree_i['mphi_phi0'], file_tree_i['A6_phi0'], file_tree_i['lambda6_phi0']), phi)))
    eps1_3.append(float(eps1_(lambda phi : V_MSSM(phi, infl_type, file_tree_i['mphi_gut'], file_tree_i['A6_gut'], file_tree_i['lambda6_gut']),phi)))

plt.figure()
plt.plot(phi_list, eps1_2, label = r'$V_{\mathrm{tree}}(p_{\mathrm{tree}})$', ls = '-', color = color_llebp1, lw = 2)
plt.plot(phi_list, eps1_1, label = r'$V_{\mathrm{llebp1}}(p_{\mathrm{llebp1}})$', ls = ':', color = color_llebp2, lw = 2)
plt.plot(phi_list, eps1_3, label = r'$V_{\mathrm{llebp1}}(p_{\mathrm{tree}})$', ls = '--', color = color_uddbp1, lw = 2)
plt.axvline(file_tree_i['phi0'], ls=':', lw=4, color='gray')
plt.axhspan(0, 1, alpha=0.25, color='orange', label='slow-roll')
plt.legend(fontsize = 25)
plt.locator_params(axis='y', nbins=6)
plt.locator_params(axis='x', nbins=6)
plt.text(float(file_tree_i['phi0'])-0.4e14, 1.5e-6, r'$\phi_0$', color = 'gray', fontsize=29)
plt.xlabel(r'$\phi$ (GeV)')
plt.ylabel(r'$\epsilon_1(\phi)$')
plt.semilogy()
plt.xlim(0,float(file_tree_i['phi0'])*1.3)
plt.grid()
define_plot_resolution()

plt.savefig('save_plots/eps1.pdf', bbox_inches='tight')
print('eps1.pdf saved')

phi_list = np.linspace(float(file_tree_i['phi0'])*0.95,float(file_tree_i['phi0'])*1.05, 200)
V1, V2, V3 = [], [], []
for phi in phi_list:
    V1.append(float(V_MSSM(phi, infl_type, file_i['mphi_gut'], file_i['A6_gut'], file_i['lambda6_gut'])))
    V2.append(float(V_MSSM(phi, 'tree', file_tree_i['mphi_phi0'], file_tree_i['A6_phi0'], file_tree_i['lambda6_phi0'])))
    V3.append(float(V_MSSM(phi, infl_type, file_tree_i['mphi_gut'], file_tree_i['A6_gut'], file_tree_i['lambda6_gut'])))

plt.figure()
plt.plot(phi_list, V2, label = r'$V_{\mathrm{tree}}(p_{\mathrm{tree}})$', ls = '-', color = color_llebp1, lw = 2)
plt.plot(phi_list, V3, label = r'$V_{\mathrm{llebp1}}(p_{\mathrm{tree}})$', ls = '--', color = color_uddbp1, lw = 2)
plt.axvline(file_tree_i['phi0'], ls=':', lw=4, color='gray')

plt.legend(fontsize = 25)
plt.locator_params(axis='y', nbins=6)
plt.locator_params(axis='x', nbins=6)
plt.xlabel(r'$\phi$ (GeV)')
plt.ylabel(r'$V(\phi)$ (GeV$^4$)')
#plt.xlim(0,float(file_tree_i['phi0'])*1.3)
#plt.ylim(bottom = 0)
plt.grid()
define_plot_resolution()

plt.savefig('save_plots/potential_zoomed.pdf', bbox_inches='tight')
print('potential_zoomed.pdf saved')
