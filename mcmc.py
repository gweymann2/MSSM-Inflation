import pandas as pd
from mpmath import *
import scipy.optimize as spo
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = 16, 10
import numpy as np
import time
from slowroll_computings import ns_v2
from def_potentials import Vrge

def log_likelihood(ns):
    return -0.5 * ((ns - 0.9665) / (3 * 0.0038)) ** 2

def mcmc(mphi_start, A6_start, lambda6_start, n):
    mphi, A6, lam6 = mphi_start, A6_start, lambda6_start
    ns = ns_v2(lambda phi: Vrge(phi, 1, mphi, A6, lam6))
    A6_list, lam6_list, mphi_list, ns_list, time_step_list = [A6], [lam6], [mphi], [ns], []
    count = 0
    for i in range(n):
        start = time.process_time()
        A6next = A6 + mpmathify('5e-16') * np.random.randn()
        lam6next = lam6 + mpmathify('5e-18') * np.random.randn()
        mphinext = mphi + mpmathify('5e-17') * np.random.randn()
        ns_next = ns_v2(lambda phi: Vrge(phi, 1, mphinext, A6next, lam6next))
        p = np.exp(log_likelihood(ns_next)) / np.exp(log_likelihood(ns))
        p_or_1 = min(p, 1)
        if np.random.uniform(0, 1) < p_or_1:
            print('\n' + str(ns) + ' to ' + str(ns_next) + ' of proba ' + str(p_or_1) + ' accepted.')
            A6, lam6, mphi, ns = A6next, lam6next, mphinext, ns_next
        else:
            print('\n' + str(ns) + ' to ' + str(ns_next) + ' of proba ' + str(p) + ' refused.')
        A6_list.append(A6)
        lam6_list.append(lam6)
        mphi_list.append(mphi)
        ns_list.append(ns)
        time_step = time.process_time() - start
        time_step_list.append(time_step)
        count += 1
        print("step " + str(count) + "/" + str(n) + " took " + str(time_step) + " sec.\n")

    return mphi_list, A6_list, lam6_list, ns_list, time_step_list



A6 = mpmathify('47072.6105683423108523859567035012099')
lam6 = mpmathify('4.49761902115191390429039481224163414e-2')
mphi = mpmathify('7446.949961062429469452587821770735257098219127860459429229360941951255418616783845167127696640252099638394991071952483716844947908')

mphi_list, A6_list, lam6_list, ns_list, time_step_list = mcmc(mphi, A6, lam6, 3)
print(mphi_list,'\n\n\n', A6_list,'\n\n\n', lam6_list,'\n\n\n', ns_list)
data = pd.DataFrame([mphi_list, A6_list, lam6_list, ns_list])
data.to_csv('mcmc3D.csv')