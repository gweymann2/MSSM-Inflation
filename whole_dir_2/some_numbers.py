import mpmath as mp
from mpmath import *
import numpy as np
import pandas as pd
import time

dps_def = 15
mp.dps = dps_def
extradps_def = 25
addprec_deriv = 30

with extradps(extradps_def):
    phigut = mp.mpf('3e16')
    lnMpinGev = mp.mpf('42.334')
    Mp = mp.exp(lnMpinGev)
    pre = mp.mpf('1')/(mp.mpf('8')*mp.pi**2)
    b1, b2, b3 = mp.mpf('33')/mp.mpf('5'), mp.mpf('1'), mp.mpf('-3')
    g1gut = mp.sqrt(mp.mpf('5')/mp.mpf('3'))*mp.mpf('5.45185741e-01')
    g2gut = mp.mpf('6.90473022e-01')
    g3gut = mp.mpf('6.84972506e-01')
    m1gut = mp.mpf('1.36108022e+02')
    m2gut = mp.mpf('1.14286222e+03')
    m3gut = mp.mpf('8.98639714e+02')
    
def say_interval(a, b, c, run=False):
    
    #if run == 'mphi_lle': fun = mphi_lle
    #elif run == 'A6_lle': fun = A6_lle
    #elif run == 'At_lle': fun = lambda phi, A6 : A6_lle(phi, (3-np.sqrt(3))/(6-np.sqrt(3))*A6)
    #elif run == 'lambda6_lle': fun = lambda6_lle
    #elif run == 'mphi_udd': fun = mphi_udd
    #elif run == 'A6_udd': fun = A6_lle
    #elif run == 'At_udd': fun = lambda phi, A6 : A6_udd(phi, (3-np.sqrt(3))/(6-np.sqrt(3))*A6)
    #elif run == 'lambda6_udd': fun = lambda6_udd
    #else: fun = lambda x, y : y
        
    #scale_end = 2000
    
    print(''+nstr(mp.mpf(str(b)), 5)+ '\ _{'+
          nstr((mp.mpf(str(a))-mp.mpf(str(b))), 5)+'}^{'+
          nstr((mp.mpf(str(c))-mp.mpf(str(b))), 5)+'}', end='\n')

    
with extradps(16):
    print('ns')
    say_interval(0.953, 0.9653, 0.9776)
    print('\n')    
"""
with extradps(16):
    print('phi0 fixed\nlle @ phi0')
    say_interval(10650.455692643742031518505129009, 11004.508211867244205999371150057, 11367.618488211001778954605382912)
    print('')
    say_interval(mp.mpf('47701.254377099194642532559017575'), mp.mpf('49287.699289601511129391481068659'), mp.mpf('50914.721470906776995961678787858'))
    print(end = ' & ')
    say_interval(0.023173525815093188847712295923457, 0.023943629521460824176232576227186, 0.024733432746701719987064099360427)
    print('\n')
"""

with extradps(16):
    print('phi0 fixed\nlle @ ewsb')
    say_interval(10665.69401807286576248665310689781181586683634445208207022797188100735931891794, 11019.25693553070001377528834128628711847218621553564936822465536360369123314877,11381.89670132814036764944415823079614699834040600792785349829948733200619199973)
    say_interval(46737.61884158651466749087001029901003911914405860648742924256339500345704265127, 48324.06375408883115434979206138301003911914405860648742924256339500345704265127, 49951.08593539409702091998978058201003911914405860648742924256339500345704265127)
    say_interval(0.043718092834228378112049517992069706609593226142907098925348223398436256803455, 0.045170934563864341339734867292430910805111555725072590594475696118626654884889, 0.046660940486889033942200384690518252529752170115853766752566842004745718229050)

with extradps(16):
    print('phi0 fixed\ntreelle @ ewsb')
    say_interval(10542.83324814976278374929209741327516117428369555805690394468506021237438246626, 10892.21467718461251647474081234452132510397534318455184070438715239623099398472,11250.56303818475744775860145602725865657981146398997300458884635760645704148137)
    say_interval(46116.40492158638914057729352835901003911914405860648742924256339500345704265127, 47681.10069518099027352913853091001003911914405860648742924256339500345704265127, 49285.81147636621068790858344124801003911914405860648742924256339500345704265127)
    say_interval(0.043077357217049264435302710082067640009844917676077328043943210643819164176533, 0.044509024609870681758890427940618720985078884923888376314523595415989076015007, 0.045977304986863061630908430698328033796986796781485740636615742609786284297514)


from scipy.interpolate import interp1d

onto_interp = [40000]

def f(_file1,_file2,_file3):

    file = pd.read_csv(_file1, dtype=str)
    
    a1 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['phi0']]), kind='cubic')(onto_interp)[0]
    b1 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['mphi_2']]), kind='cubic')(onto_interp)[0]
    c1 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['A6_2']]), kind='cubic')(onto_interp)[0]
    d1 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['lambda6_2']]), kind='cubic')(onto_interp)[0]

    file = pd.read_csv(_file2, dtype=str)
    a2 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['phi0']]), kind='cubic')(onto_interp)[0]
    b2 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['mphi_2']]), kind='cubic')(onto_interp)[0]
    c2 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['A6_2']]), kind='cubic')(onto_interp)[0]
    d2 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['lambda6_2']]), kind='cubic')(onto_interp)[0]

    file = pd.read_csv(_file3, dtype=str)
    a3 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['phi0']]), kind='cubic')(onto_interp)[0]
    b3 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['mphi_2']]), kind='cubic')(onto_interp)[0]
    c3 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['A6_2']]), kind='cubic')(onto_interp)[0]
    d3 = interp1d(np.array([float(x) for x in file['A6_gut']]), np.array([float(x) for x in file['lambda6_2']]), kind='cubic')(onto_interp)[0]

    say_interval(a1,a2,a3)   
    say_interval(b1,b2,b3)   
    say_interval(c1,c2,c3)   
    say_interval(d1,d2,d3)   
    
list_files = [
'save_file/lle_bp1_0_3047_9627.csv', 
'save_file/lle_bp1_0_3047_9665.csv', 
'save_file/lle_bp1_0_3047_9703.csv', 
'save_file/tree_0_3047_9627_high_llebp1.csv', 
'save_file/tree_0_3047_9665_high_llebp1.csv',
'save_file/tree_0_3047_9703_high_llebp1.csv', 
 ]

print('')
f(list_files[0],list_files[1],list_files[2],)
print('')
f(list_files[3],list_files[4],list_files[5],)

