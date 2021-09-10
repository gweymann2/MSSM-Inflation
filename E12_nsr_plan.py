from E0_aspic import x_rrad, norm_eps1, norm_eps2
import mpmath as mp
import matplotlib.pyplot as plt
from A_def_potentials import Mp
from getdist.mcsamples import loadMCSamples
import getdist.plots as gdplt
import numpy as np
import pandas as pd

Pstar = mp.mpf('2.105e-9')

def ns_r_plan(phi0, alpha, lnRrad):
    xstar = x_rrad(alpha, phi0/Mp, lnRrad, Pstar)
    eps1 = norm_eps1(xstar,alpha,phi0/Mp)
    eps2 = norm_eps2(xstar,alpha,phi0/Mp)
    ns = 1-2*eps1-eps2
    r = 16*eps1
    return float(ns), float(r)

def generate_examples(phi0_start, alpha_start, delta_alpha):
    L = []
    for i in range(-18,32):
        print('.',end="", flush=True)
        alpha = alpha_start+i*delta_alpha
        L.append([*ns_r_plan(phi0_start, alpha, 0), float(mp.log(mp.fabs(1-alpha)))])
    return L

delta_alpha = [0.00000000000006*(0.0000000006/0.00000000000006)**(i/50) for i in range(50)]


phi0_start = np.array([mp.mpf(phi0) for phi0 in np.array(pd.read_csv('E11_tree09665Inst.csv', dtype = str, index_col=0)['phi0B'])])
alpha_start = np.array([mp.mpf(phi0) for phi0 in np.array(pd.read_csv('E11_tree09665Inst.csv', dtype = str, index_col=0)['alpha'])])

example_tot = []
for i in range(50):
    example_i = generate_examples(phi0_start[i], alpha_start[i], delta_alpha[i])
    example_tot = example_tot+example_i
    print(' '+str(i+1)+'/50')
example_tot = np.transpose(example_tot)


s = loadMCSamples('plikHM_TTTEEE_lowl_lowE_BK15_lensing/base_r_plikHM_TTTEEE_lowl_lowE_BK15_lensing')
plot = gdplt.get_subplot_plotter(subplot_size=4)
plot.triangle_plot(s, ['A','ns', 'r'], shaded = True)
nsr = plot.subplots[2,1]

nsr.semilogy()
nsr.axis(ymin=10**-20,ymax=100)
nsr.scatter(example_tot[0], example_tot[1], c=example_tot[2], cmap = 'gist_rainbow', marker = 's', s=1000)
plot.fig.savefig('E13_triplot.png')
