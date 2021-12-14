from E_contour import *

mp.dps = 500
mp.prec = 166

lnMpinGev = mp.mpf('42.334')
Mp = mp.exp(lnMpinGev)

def define_plot_resolution():
    fig = plt.gcf()  # get current figure

    DPI = fig.get_dpi()
    #     fig.set_size_inches(1920.0 / float(DPI), 1080.0 / float(DPI))
    fig.set_size_inches(12, 8)
    ax = plt.gca()
    for tickLabel in ax.get_xticklabels() + ax.get_yticklabels():
        tickLabel.set_fontsize(29)
    ax.yaxis.label.set_size(29)
    ax.xaxis.label.set_size(29)
    ax.yaxis.offsetText.set_fontsize(29)
    ax.xaxis.offsetText.set_fontsize(29)
    return


m0_Rge_09530 = pd.read_csv('m0_LLe_9530.csv', engine='python', dtype=str)
m0_Rge_09653 = pd.read_csv('m0_lle_9653.csv', engine='python', dtype=str)
m0_Rge_09776 = pd.read_csv('m0_lle_9776.csv', engine='python', dtype=str)
m0_Tree_09530 = pd.read_csv('20pts_m0_Tree_lle_9530.csv', engine='python', dtype=str)
m0_Tree_09653 = pd.read_csv('20pts_m0_Tree_lle_9653.csv', engine='python', dtype=str)
m0_Tree_09776 = pd.read_csv('20pts_m0_Tree_lle_9776.csv', engine='python', dtype=str)


# m10_Rge_09530 = pd.read_csv('20pts_m10_Rge_09530.csv',engine='python',dtype=str)
# m10_Rge_09653 = pd.read_csv('20pts_m10_Rge_09653.csv',engine='python',dtype=str)
# m10_Rge_09776 = pd.read_csv('20pts_m10_Rge_09776.csv',engine='python',dtype=str)
# m10_Tree_09530 = pd.read_csv('20pts_m10_Tree_09530.csv',engine='python',dtype=str)
# m10_Tree_09653 = pd.read_csv('20pts_m10_Tree_09653.csv',engine='python',dtype=str)
# m10_Tree_09776 = pd.read_csv('20pts_m10_Tree_09776.csv',engine='python',dtype=str)
# print(m10_Tree_09776['phi0B'])

def pd_to_array(file, column, dtype=float):
    if dtype == float:
        return np.array([mp.mpf(x) for x in file[column]], dtype=float)
    else:
        return np.array([mp.mpf(x) for x in file[column]])


print(end='| ')
for col in m0_Rge_09530.columns[1:]:
    print(col, end=' | ')

plt.figure(1)
plt.plot(pd_to_array(m0_Tree_09653, 'phi0B'), pd_to_array(m0_Tree_09653, 'mphi'),'blue')
plt.plot(pd_to_array(m0_Rge_09653, 'phi0B'), pd_to_array(m0_Rge_09653, 'mphi'),'red')
plt.fill_between(pd_to_array(m0_Tree_09530, 'phi0B'), pd_to_array(m0_Tree_09530, 'mphi'), pd_to_array(m0_Tree_09776, 'mphi'), color = 'blue', alpha=0.2, label='Tree')
plt.fill_between(pd_to_array(m0_Rge_09530, 'phi0B'), pd_to_array(m0_Rge_09530, 'mphi'), pd_to_array(m0_Rge_09776, 'mphi'), color='red', alpha=0.2, label='Rge')
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'$m_{\phi,GUT}$ (GeV)')
plt.legend(fontsize=25)
define_plot_resolution()

plt.figure(2)
plt.plot(pd_to_array(m0_Tree_09653, 'phi0B'), pd_to_array(m0_Tree_09653, 'A6'),'blue')
plt.plot(pd_to_array(m0_Rge_09653, 'phi0B'), pd_to_array(m0_Rge_09653, 'A6'),'red')
plt.fill_between(pd_to_array(m0_Tree_09530, 'phi0B'), pd_to_array(m0_Tree_09530, 'A6'), pd_to_array(m0_Tree_09776, 'A6'), color = 'blue', alpha=0.2, label='Tree')
plt.fill_between(pd_to_array(m0_Rge_09530, 'phi0B'), pd_to_array(m0_Rge_09530, 'A6'), pd_to_array(m0_Rge_09776, 'A6'), color='red', alpha=0.2, label='Rge')
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'$A_{6,GUT}$ (GeV)')
plt.legend(fontsize=25)
define_plot_resolution()

plt.figure(3)
plt.plot(pd_to_array(m0_Tree_09653, 'phi0B'), pd_to_array(m0_Tree_09653, 'lambda6'),'blue')
plt.plot(pd_to_array(m0_Rge_09653, 'phi0B'), pd_to_array(m0_Rge_09653, 'lambda6'),'red')
plt.fill_between(pd_to_array(m0_Tree_09530, 'phi0B'), pd_to_array(m0_Tree_09530, 'lambda6'), pd_to_array(m0_Tree_09776, 'lambda6'), color = 'blue', alpha=0.2, label='Tree')
plt.fill_between(pd_to_array(m0_Rge_09530, 'phi0B'), pd_to_array(m0_Rge_09530, 'lambda6'), pd_to_array(m0_Rge_09776, 'lambda6'), color='red', alpha=0.2, label='Rge')
plt.xlabel(r'$\phi_0$ (GeV)')
plt.ylabel(r'$\lambda_{6,GUT}$')
plt.legend(fontsize=25)
plt.semilogy()
define_plot_resolution()

plt.figure(4)
plt.plot(pd_to_array(m0_Tree_09653, 'mphi'), pd_to_array(m0_Tree_09653, 'A6'),'blue')
plt.plot(pd_to_array(m0_Rge_09653, 'mphi'), pd_to_array(m0_Rge_09653, 'A6'),'red')
plt.fill_between([], [], [], color = 'blue', alpha=0.2, label='Tree')
plt.fill_between([], [], [], color='red', alpha=0.2, label='Rge')
plt.xlabel(r'$A_6$ (GeV)')
plt.ylabel(r'$m_{\phi,GUT}$ (GeV)')
plt.legend(fontsize=25)
define_plot_resolution()

