from mpmath import *
import pandas as pd
from D0_mcmc import mcmc
import corner
import numpy as np
import matplotlib.pyplot as plt

### EXECUTE MCMC ###
nbr_it = 5
A6_start, lam6_start = mpmathify('47072.6105683423108523859567035012099'), mpmathify('4.49761902115191390429039481224163414e-2')
A6_list, lam6_list, ns_list, r_list, eps1_list, eps2_list, eps3_list, alpha_list, time_step_list = mcmc(nbr_it, A6_start, lam6_start)

### WRITE CHAIN IN A FILE ###
data = pd.DataFrame(data={'A6':A6_list, 'lam6':lam6_list, 'ns':ns_list, 'r':r_list, 'eps1':eps1_list, 'eps2':eps2_list, 'eps3':eps3_list, 'alpha':alpha_list})
name_file = 'D2_results_mcmc.csv'
data.to_csv(name_file)
#df = pd.read_csv(name_file, dtype = str, float_precision=5, index_col=0)
#print(df)

### CORNER ###
A6_list_float = []
for el in A6_list:
    A6_list_float.append(float(el-A6_start))
lam6_list_float = []
for el in lam6_list:
    lam6_list_float.append(float(el-lam6_start))
my_samples = np.transpose(np.array([np.array(lam6_list_float)*1.25e17, np.array(A6_list_float)*1.25e15, ns_list]))
plt.figure(1)
figure = corner.corner(my_samples, labels = [r'$(\lambda_6-\lambda_{6centre})*1.25e17$', r'$(A_6-A_{6centre})*1.25e15$', "ns"])
plt.savefig('D3_corner')
