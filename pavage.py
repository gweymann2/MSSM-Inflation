from mpmath import *
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = 16, 10
import numpy as np
import time
from slowroll_computings import ns_v2
from def_potentials import Vrge

start_tot = time.process_time()

mphi = mpmathify('7446.949961062429469452587821770735257098219127860459429229360941951255418616783845167127696640252099638394991071952483716844947908')
A6 = mpmathify('47072.6105683423108523859567035012099')
lambda6 = mpmathify('4.49761902115191390429039481224163414e-2')
count = 0
side = 0 # vrai side est 2 side +1
step_A = mpmathify('0.8e-15')
step_lambda = mpmathify('0.8e-17')
grid = np.zeros((2*side+1,2*side+1))
for j in range(-side, side+1):
    A6_j = A6+j*step_A
    for k in range(-side, side+1):
        start = time.process_time()
        lambda6_k = lambda6+k*step_lambda
        grid[j+side][k+side] = ns_v2(lambda phi : Vrge(phi, 1, mphi, A6_j, lambda6_k))
        print('\n',grid)
        count += 1
        print("step "+ str(count) +"/"+str((2*side+1)**2)+" took "+str(time.process_time() - start) +" sec ",'\n')
np.savetxt('[A lambda] grid', grid, delimiter = ', ')
print("total time for "+str((2*side+1)**2)+" steps took "+str(time.process_time() - start_tot) +" sec ")

plt.figure(1)
plt.imshow(grid, extent=[-8.5,8.5,8.5,-8.5], cmap='YlOrRd')
cb = plt.colorbar()
cb.set_label(label='ns', fontsize=18)
plt.xlabel(r'$(\lambda_6-\lambda_{6centre})*1.25e17$', fontsize = 18)
plt.ylabel(r'$(A_6-A_{6centre})*1.25e15$', fontsize = 18)
plt.savefig('[A lambda] ns distrib')

plt.figure(2)
palette = plt.matplotlib.colors.LinearSegmentedColormap('jet3',plt.cm.datad['jet'],2048)
palette.set_under(alpha=0.0)
plt.imshow(np.log(((np.array(grid)-0.9665)/0.0038)**2), cmap = palette)
cb = plt.colorbar()
cb.set_label(label='log((ns-0.9665/0.0038)²)', fontsize=18)
plt.xlabel(r'$(\lambda_6-\lambda_{6centre})*1.25e17+8$', fontsize = 18)
plt.ylabel(r'$(A_6-A_{6centre})*1.25e15+8$', fontsize = 18)
plt.savefig('[A lambda] chi2 distrib')
