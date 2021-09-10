import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def define_plot_resolution():

    fig = plt.gcf()

    DPI = fig.get_dpi()
    fig.set_size_inches(1920.0 / float(DPI), 1080.0 / float(DPI))
    ax = plt.gca()
    for tickLabel in ax.get_xticklabels()+ax.get_yticklabels():
        tickLabel.set_fontsize(35)
    ax.yaxis.label.set_size(35)
    ax.xaxis.label.set_size(35)
    ax.yaxis.offsetText.set_fontsize(35)
    ax.xaxis.offsetText.set_fontsize(35)
    return

def pol_non_entier(x, A, u, B, v, C, w, D, t):
    return A*x**u+B*x**v+C*x**w+D*x**t

def sci(a):
    return("%.5e"%a)

aspic09551lnRradm10 = pd.read_csv('aspic_data/gmssmi_september_09551_lnRradm10.dat',delimiter="  ",names = ['phi0B','mphi','r','ns','Treh','lambdaB','AB'],engine='python')#,'phi0B','lambdaB','AB'
aspic09779lnRradm10 = pd.read_csv('aspic_data/gmssmi_september_09779_lnRradm10.dat',delimiter="  ",names = ['phi0B','mphi','r','ns','Treh','lambdaB','AB'],engine='python')#,'phi0B','lambdaB','AB'
aspic09779lnRrad0 = pd.read_csv('aspic_data/gmssmi_september_09779_lnRrad0.dat',delimiter="  ",names = ['phi0B','mphi','r','ns','Treh','lambdaB','AB'],engine='python')#,'phi0B','lambdaB','AB'
aspic09551lnRrad0 = pd.read_csv('aspic_data/gmssmi_september_09551_lnRrad0.dat',delimiter="  ",names = ['phi0B','mphi','r','ns','Treh','lambdaB','AB'],engine='python')#,'phi0B','lambdaB','AB'

python09551lnRradm10 = pd.read_csv('E5_tree09551m10.csv', dtype = float, float_precision=5, index_col=0)
python09779lnRradm10 = pd.read_csv('E4_tree09779m10.csv', dtype = float, float_precision=5, index_col=0)
python09779lnRrad0 = pd.read_csv('E2_tree09779Inst.csv', dtype = float, float_precision=5, index_col=0)
python09551lnRrad0 = pd.read_csv('E3_tree09551Inst.csv', dtype = float, float_precision=5, index_col=0)

plt.figure(0)
plt.fill_between(python09551lnRradm10['phi0B'],python09551lnRradm10['mphi'],python09779lnRradm10['mphi'], label='0.9551 < ns < 0.9779 & lnRrad = -10')
plt.fill_between(python09551lnRrad0['phi0B'],python09551lnRrad0['mphi'],python09779lnRrad0['mphi'], label='0.9551 < ns < 0.9779 & lnRrad = 0')

x = np.array(aspic09551lnRradm10['phi0B'])
linestyles = ['-','-.',':','--']
for i, name_file in enumerate([python09551lnRradm10, python09779lnRradm10, python09551lnRrad0, python09779lnRrad0]) :
    power = []
    n = 4
    y = np.array(name_file['mphi'])
    erry = y/10
    s = len(x)
    for k in range (n):
        xs,ys=x[k*s//n:(k+1)*s//n],y[k*s//n:(k+1)*s//n]
        power.append((np.log(ys[-1])-np.log(ys[0]))/(np.log(xs[-1])-np.log(xs[0])))
    par, par_va = curve_fit(lambda x, A, B, C, D : pol_non_entier(x, A, power[0], B, power[1], C, power[2], D, power[3]), x, y, p0=[1e-25,1e-25,1e-25,1e-25], sigma=erry, absolute_sigma=True)
    plt.plot(x, pol_non_entier(x,par[0],power[0],par[1],power[1],par[2],power[2],par[3],power[3]), linestyle=linestyles[i], color='black', label=r'$m_{\phi}$ = '+str(sci(par[0]))+' * '+r'$\phi_0$'+' ^ '+str(round(power[0],4))+' + '+str(sci(par[1]))+' * '+r'$\phi_0$'+' ^ '+str(round(power[1],4))+' + '+str(sci(par[2]))+' * '+r'$\phi_0$'+' ^ '+str(round(power[2],4))+' + '+str(sci(par[3]))+' * '+r'$\phi_0$'+' ^ '+str(round(power[3],4)))#+' * (1+'+str(sci(c0/1.2))+' * '+r'$\phi_0$'+' ^ '+str(sci(d0))+') * (1+'+str(sci(e0/0.15))+' * '+r'$\phi_0$'+' ^ '+str(sci(f0))+')')

plt.legend(fontsize=22)
plt.xlabel(r'$\phi_0$')
plt.ylabel(r'$m\phi_{python}$')
plt.loglog()
define_plot_resolution()
plt.savefig('E7_phi0_mphi_loglog')


plt.figure(1)
plt.plot(aspic09551lnRradm10['phi0B'],(aspic09551lnRradm10['mphi']-python09551lnRradm10['mphi'])/python09551lnRradm10['mphi'], label='ns = 0.9551 & lnRrad = -10')
plt.plot(aspic09779lnRradm10['phi0B'],(aspic09779lnRradm10['mphi']-python09779lnRradm10['mphi'])/python09779lnRradm10['mphi'], label='ns = 0.9779 & lnRrad = -10')
plt.plot(aspic09551lnRrad0['phi0B'],(aspic09551lnRrad0['mphi']-python09551lnRrad0['mphi'])/python09551lnRrad0['mphi'], label='ns = 0.9551 & lnRrad = 0')
plt.plot(aspic09779lnRrad0['phi0B'],(aspic09779lnRrad0['mphi']-python09779lnRrad0['mphi'])/python09779lnRrad0['mphi'], label='ns = 0.9779 & lnRrad = 0')
plt.legend(fontsize=25)
plt.xlabel(r'$\phi_0$')
plt.ylabel(r'$\frac{m\phi_{aspic}-m\phi_{python}}{m\phi_{python}}$')
plt.semilogx()
define_plot_resolution()
plt.savefig('E8_relative_error_python_aspic')

plt.figure(2)
plt.plot(python09551lnRradm10['phi0B'],(python09551lnRradm10['ns']-0.9551)/0.9551, label='ns = 0.9551 & lnRrad = -10')
plt.plot(python09779lnRradm10['phi0B'],(python09779lnRradm10['ns']-0.9779)/0.9779, label='ns = 0.9779 & lnRrad = -10')
plt.plot(python09551lnRrad0['phi0B'],(python09551lnRrad0['ns']-0.9551)/0.9551, label='ns = 0.9551 & lnRrad = 0')
plt.plot(python09779lnRrad0['phi0B'],(python09779lnRrad0['ns']-0.9779)/0.9779, label='ns = 0.9779 & lnRrad = 0')
plt.legend(fontsize=25)
plt.title('python',fontsize=25)
plt.xlabel(r'$\phi_0$')
plt.ylabel(r'$\frac{ns_{out}-ns_{in}}{ns_{in}}$')
plt.semilogx()
define_plot_resolution()
plt.savefig('E9_consistency_check_python')


plt.figure(3)
plt.plot(aspic09551lnRradm10['phi0B'],(aspic09551lnRradm10['ns']-0.9551)/0.9551, label='ns = 0.9551 & lnRrad = -10')
plt.plot(aspic09779lnRradm10['phi0B'],(aspic09779lnRradm10['ns']-0.9779)/0.9779, label='ns = 0.9779 & lnRrad = -10')
plt.plot(aspic09551lnRrad0['phi0B'],(aspic09551lnRrad0['ns']-0.9551)/0.9551, label='ns = 0.9551 & lnRrad = 0')
plt.plot(aspic09779lnRrad0['phi0B'],(aspic09779lnRrad0['ns']-0.9779)/0.9779, label='ns = 0.9779 & lnRrad = 0')
plt.legend(fontsize=25)
plt.title('aspic',fontsize=25)
plt.xlabel(r'$\phi_0$')
plt.ylabel(r'$\frac{ns_{out}-ns_{in}}{ns_{in}}$')
plt.semilogx()
define_plot_resolution()
plt.savefig('E10_consistency_check_aspic')

plt.figure(4)
plt.plot(aspic09551lnRradm10['phi0B'],0.9551*(aspic09551lnRradm10['mphi']-python09551lnRradm10['mphi'])/(python09551lnRradm10['mphi']*(aspic09551lnRradm10['ns']-0.9551)), label='ns = 0.9551 & lnRrad = -10')
plt.plot(aspic09779lnRradm10['phi0B'],0.9779*(aspic09779lnRradm10['mphi']-python09779lnRradm10['mphi'])/(python09779lnRradm10['mphi']*(aspic09779lnRradm10['ns']-0.9779)), label='ns = 0.9779 & lnRrad = 0')
plt.plot(aspic09551lnRrad0['phi0B'],0.9551*(aspic09551lnRrad0['mphi']-python09551lnRrad0['mphi'])/(python09551lnRrad0['mphi']*(aspic09551lnRrad0['ns']-0.9551)), label='ns = 0.9551 & lnRrad = 0')
plt.plot(aspic09779lnRrad0['phi0B'],0.9779*(aspic09779lnRrad0['mphi']-python09779lnRrad0['mphi'])/(python09779lnRrad0['mphi']*(aspic09779lnRrad0['ns']-0.9779)), label='ns = 0.9779 & lnRrad = 0')
plt.legend(fontsize=25)
plt.xlabel(r'$\phi_0$')
plt.ylabel(r'$\frac{m\phi_{aspic}-m\phi_{python}}{m\phi_{python}}*\frac{ns_{in,aspic}}{ns_{out,aspic}-ns_{in,aspic}}$')
plt.semilogx()
define_plot_resolution()
plt.savefig('E11_ok')


