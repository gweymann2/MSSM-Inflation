import mpmath as mp
from A_def_potentials import Vrge
from B0_slowroll_computings import N, eps1f, eps2f, eps3f, endinf
from B0_slowroll_computings import Pstar
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = 16, 10

mphi_ex = mp.mpf('7446.949961062429469357133726782183707876898317924005895581087450892568320619984020659573100501291276399015423639649744434796823819')
A6_ex = mp.mpf('47072.6105683423108523859567035012099')
lambda6_ex = mp.mpf('4.49761902115191390429039481224163414e-2')

V = lambda phi : Vrge(phi, 1, mphi_ex, A6_ex, lambda6_ex)

phi_end = endinf(V, 9.3e14)

kstar, lnMpcToKappa, HubbleSquareRootOf3OmegaRad, RelatDofRatio = 0.05, 130.282, 7.5437e-63, 1
N0 = mp.log(kstar) - lnMpcToKappa - 0.5*mp.log(HubbleSquareRootOf3OmegaRad) -0.25*mp.log(RelatDofRatio)
lnRrad = 0
la_func = lambda phi : lnRrad - N0 - 0.25*mp.log(9/eps1f(V, phi)/(3-eps1f(V, phi_end))*V(phi_end)/V(phi))+0.25*mp.log(8*mp.pi**2*Pstar)

mp.plot([lambda phi : N(V,phi, phi_end), la_func],  xlim=[9.131718185e14,9.131718208e14], ylim=None, points=40, file=None, dpi=None, singularities=[], axes=None)
mp.plot([lambda phi : mp.log(eps1f(V,phi))],  xlim=[9.131718185e14,9.131718208e14], ylim=None, points=1000, file=None, dpi=None, singularities=[], axes=None)
mp.plot([lambda phi : eps2f(V,phi)],  xlim=[9.131718185e14,9.131718208e14], ylim=None, points=1000, file=None, dpi=None, singularities=[], axes=None)
mp.plot([lambda phi : eps3f(V,phi)],  xlim=[9.131718185e14,9.131718208e14], ylim=None, points=1000, file=None, dpi=None, singularities=[], axes=None)
