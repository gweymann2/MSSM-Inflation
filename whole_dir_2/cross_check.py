import mpmath as mp
from mpmath import *
import numpy as np
import pandas as pd
import time
print('module correctly imported')

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

def DV(V,phi):
    return mp.diff(V, phi, addprec=addprec_deriv)

def D2V(V, phi):
    return mp.diff(V, phi, 2, addprec=addprec_deriv)
def eps1_(V, phi):
    return Mp**2/2*(DV(V,phi)/V(phi))**2

def eps2_(V, phi):
    return 2*Mp**2*((DV(V,phi)/V(phi))**2-D2V(V,phi)/V(phi))


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

dps_def = 32
mp.dps = dps_def
extradps_def = 25
addprec_deriv = 30
mp.dps = dps_def
with extradps(extradps_def):
    phigut = mp.mpf('3e16')
    lnMpinGev = mp.mpf('42.334')
    Mp = mp.exp(lnMpinGev)
    pre = mp.mpf('1')/(mp.mpf('8')*mp.pi**2)
    b1, b2, b3 = mp.mpf('33')/mp.mpf('5'), mp.mpf('1'), mp.mpf('-3')
    bp = ''
    fac_corr = mp.sqrt(mp.mpf('2'))
    fac_highscale = mp.mpf('1')
    
g1_ = lambda phi, g1gut : g1gut/(mp.sqrt(1-pre*b1*g1gut**2*mp.log(phi/phigut)))
g2_ = lambda phi, g2gut : g2gut/(mp.sqrt(1-pre*b2*g2gut**2*mp.log(phi/phigut)))
g3_ = lambda phi, g3gut : g3gut/(mp.sqrt(1-pre*b3*g3gut**2*mp.log(phi/phigut)))
M1_ = lambda phi, m1gut, g1gut : m1gut*(g1_(phi, g1gut)/g1gut)**mp.mpf('2')
M2_ = lambda phi, m2gut, g2gut : m2gut*(g2_(phi, g2gut)/g2gut)**mp.mpf('2')
M3_ = lambda phi, m3gut, g3gut : m3gut*(g3_(phi, g3gut)/g3gut)**mp.mpf('2')
mphi_run_lle_ = lambda phi_start, phi_end, mphi_start, m1gut, g1gut, m2gut, g2gut : mp.sqrt(mphi_start**2+(M2_(phi_start, m2gut, g2gut)**2-M2_(phi_end, m2gut, g2gut)**2)+mp.mpf('1')/11*(M1_(phi_start, m1gut, g1gut)**2-M1_(phi_end, m1gut, g1gut)**2))
A6_run_lle_ = lambda phi_start, phi_end, A6_start, m1gut, g1gut, m2gut, g2gut : A6_start-mp.mpf('6')*(M2_(phi_start, m2gut, g2gut)-M2_(phi_end, m2gut, g2gut))-mp.mpf('6')/11*(M1_(phi_start, m1gut, g1gut)-M1_(phi_end, m1gut, g1gut))
lambda6_run_lle_ = lambda phi_start, phi_end, lambda6_start, g1gut, g2gut : lambda6_start*(g2_(phi_start, g2gut)/g2_(phi_end, g2gut))**mp.mpf('6')*(g1_(phi_start, g1gut)/g1_(phi_end, g1gut))**(mp.mpf('6')/11)
mphi_run_udd_ = lambda phi_start, phi_end, mphi_start, m1gut, g1gut, m3gut, g3gut : mp.sqrt(mphi_start**2-mp.mpf('8')/9*(M3_(phi_start, m3gut, g3gut)**2-M3_(phi_end, m3gut, g3gut)**2)+mp.mpf('4')/99*(M1_(phi_start, m1gut, g1gut)**2-M1_(phi_end, m1gut, g1gut)**2))
A6_run_udd_ = lambda phi_start, phi_end, A6_start, m1gut, g1gut, m3gut, g3gut : A6_start+mp.mpf('16')/3*(M3_(phi_start, m3gut, g3gut)-M3_(phi_end, m3gut, g3gut))-mp.mpf('8')/33*(M1_(phi_start, m1gut, g1gut)-M1_(phi_end, m1gut, g1gut))
lambda6_run_udd_ = lambda phi_start, phi_end, lambda6_start, g1gut, g3gut : lambda6_start*(g3_(phi_start, g3gut)/g3_(phi_end, g3gut))**(mp.mpf('6')*(-mp.mpf('8'))/9)*(g1_(phi_start, g1gut)/g1_(phi_end, g1gut))**(4*mp.mpf('6')/99)
mphi_run_tree_ = lambda phi_start, phi_end, mphi_start : mphi_start
A6_run_tree_ = lambda phi_start, phi_end, A6_start : A6_start
lambda6_run_tree_ = lambda phi_start, phi_end, lambda6_start : lambda6_start

def V_MSSM_run(phi, phi_start, infl_type, mphigut, A6gut, lambda6gut):


    if infl_type == 0 or infl_type == 'tree':
        mphi_func, A6_func, lambda6_func = mphi_tree_, A6_tree_, lambda6_tree_
    elif infl_type == 1 or infl_type == 'lle':
        mphi_func = lambda phi_start, phi_end, mphi_start : mphi_run_lle_(phi_start, phi_end, mphi_start, m1gut, g1gut, m2gut, g2gut)
        A6_func = lambda phi_start, phi_end, A6_start : A6_run_lle_(phi_start, phi_end, A6_start, m1gut, g1gut, m2gut, g2gut)
        lambda6_func = lambda phi_start, phi_end, lambda6_start : lambda6_run_lle_(phi_start, phi_end, lambda6_start, g1gut, g2gut)
    elif infl_type == 2 or infl_type == 'udd':
        mphi_func = lambda phi_start, phi_end, mphi_start : mphi_run_udd_(phi_start, phi_end, mphi_start, m1gut, g1gut, m3gut, g3gut)
        A6_func = lambda phi_start, phi_end, A6_start : A6_run_udd_(phi_start, phi_end, A6_start, m1gut, g1gut, m3gut, g3gut)
        lambda6_func = lambda phi_start, phi_end, lambda6_start : lambda6_run_udd_(phi_start, phi_end, lambda6_start, g1gut, g3gut)
    else:
        return 'Error: unknown type of inflation'
    lambda6 = lambda6_func(phi_start, phi, lambda6gut)
    mphi = mphi_func(phi_start, phi,mphigut)
    A6 = A6_func(phi_start, phi,A6gut)
    V = mp.mpf('0.5')*mphi**mp.mpf('2')*phi**mp.mpf('2')-fac_corr*fac_highscale*lambda6*A6/(mp.mpf('6')*Mp**mp.mpf('3'))*phi**mp.mpf('6')+fac_highscale**2*lambda6**mp.mpf('2')*phi**mp.mpf('10')/Mp**mp.mpf('6')
    return V
mphi_run_lle = lambda phi_start, phi_end, mphi_start : mp.sqrt(mphi_start**2+(M2(phi_start)**2-M2(phi_end)**2)+mp.mpf('1')/11*(M1(phi_start)**2-M1(phi_end)**2))
lambda6_run_lle = lambda phi_start, phi_end, lambda6_start : lambda6_start*(g2(phi_start)/g2(phi_end))**mp.mpf('6')*(g1(phi_start)/g1(phi_end))**(mp.mpf('6')/11)
mphi_run_udd = lambda phi_start, phi_end, mphi_start : mp.sqrt(mphi_start**2-mp.mpf('8')/9*(M3(phi_start)**2-M3(phi_end)**2)+mp.mpf('4')/99*(M1(phi_start)**2-M1(phi_end)**2))
lambda6_run_udd = lambda phi_start, phi_end, lambda6_start : lambda6_start*(g3(phi_start)/g3(phi_end))**(mp.mpf('6')*(-mp.mpf('8'))/9)*(g1(phi_start)/g1(phi_end))**(4*mp.mpf('6')/99)
mphi_run_tree = lambda phi_start, phi_end, mphi_start : mphi_start
lambda6_run_tree = lambda phi_start, phi_end, lambda6_start : lambda6_start


def FI(phi, A6gut, dv, infl_type, at_phi0=False, potential_at_phi0=False, return_Acal = False, return_alpha=False):

    if infl_type == 'lle':
        A6 = A6_lle(phi, A6gut)
        mphi_run = mphi_run_lle
        lambda6_run = lambda6_run_lle
        ga = g1
        gb = g2
        Ma = M1
        Mb = M2
        with extradps(extradps_def):
            Ba=b1/(16*mp.pi**2)
            Bb=b2/(16*mp.pi**2)
            ea=mp.mpf('1')/mp.mpf('11')*b1/(4*mp.pi**2) 
            eb=b2/(mp.mpf('4')*mp.pi**2)
            
    elif infl_type == 'udd':
        A6 = A6_udd(phi, A6gut)
        mphi_run = mphi_run_udd
        lambda6_run = lambda6_run_udd
        ga = g1
        gb = g3
        Ma = M1
        Mb = M3
        with extradps(0):
            Ba=b1/(16*mp.pi**2)
            Bb=b3/(16*mp.pi**2)
            ea=mp.mpf('4')/mp.mpf('99')*b1/(4*mp.pi**2) 
            eb=mp.mpf('-8')/mp.mpf('9')*b3/(mp.mpf('4')*mp.pi**2)
    
    elif infl_type=='tree':
        A6 = A6_tree(phi, A6gut)
        mphi_run, lambda6_run = mphi_run_tree, lambda6_run_tree
        ga, gb, Ma, Mb = g1, g3, M1, M3
        with extradps(extradps_def):
            Ba, Bb, ea, eb = 0, 0, 0, 0
    
    else: 
        print('this infl_type is not implemented')
#     print(ea, eb, ga(phi), gb(phi), Ma(phi), Mb(phi))

    Ca=2*Ba
    Cb=2*Bb
    Fa=3*ea
    Fb=3*eb
    Da=Fa/2
    Db=Fb/2
    
    
    betaga = 2*Ba*ga(phi)**4
    betagb = 2*Bb*gb(phi)**4
    betaMa =  Ca*ga(phi)**2*Ma(phi)
    betaMb = Cb*gb(phi)**2*Mb(phi)
    betam = -ea*ga(phi)**2*Ma(phi)**2-eb*gb(phi)**2*Mb(phi)**2
    betaA = Fa*ga(phi)**2*Ma(phi)+Fb*gb(phi)**2*Mb(phi)
    betalambda = -(Da*ga(phi)**2+Db*gb(phi)**2)
    CST1p = 2*dv/phi-betam
    CST2p = -2*((Ba+Ca)*ea*Ma(phi)**2*ga(phi)**4+(Bb+Cb)*eb*Mb(phi)**2*gb(phi)**4)+3*betam #OK
    c1 = 12*(5+betalambda) #OK
    c2 = 12*(Da*betaga+Db*betagb-(5+betalambda)*(9+2*betalambda))
    c3 = fac_corr * (betaA+A6*(6+betalambda))
    c4 = fac_corr * (-(2*Ba+Ca)*Fa*Ma(phi)*ga(phi)**4-(2*Bb+Cb)*Fb*Mb(phi)*gb(phi)**4-betaA*(mp.mpf('11')+2*betalambda)+A6*(Da*betaga+Db*betagb-(5+betalambda)*(6+betalambda)))
    COEF1 = 3*(CST2p*c1 - CST1p*c2)**2 + (CST2p*c3 - CST1p*c4)*(-c2*c3 + c1*c4) #OK
    COEF2 = mp.mpf('1')/10*(6*(c1 + c2)*(CST2p*c1 - CST1p*c2) + (c3 + c4)*(-c2*c3 + c1*c4))
    COEF3 = mp.mpf('3')/100*(c1 + c2)**2 #OK
    r1 = COEF1/COEF3
    r2 = COEF2/COEF3
    Acal = mp.sqrt(-r2+mp.sqrt(-4*r1+r2**2))
    mphi_FI_phi0 = Acal/mp.sqrt(mp.mpf('40'))
    lambda6_FI_phi0 = 1/fac_highscale*(3*Mp**3/phi**4*(CST2p*c1 - CST1p*c2 + 2*mphi_FI_phi0**2*(c1 + c2))/(c2*c3 - c1*c4))
    

    if at_phi0:
        return mphi_FI_phi0, lambda6_FI_phi0, A6
    if potential_at_phi0:
        return mp.mpf('0.5')*mphi_FI_phi0**2*phi**2-fac_corr*A6*lambda6_FI_phi0*phi**6/(6*Mp**3)+lambda6_FI_phi0**2/Mp**6*phi**10
        
    with extradps(extradps_def):
        mphi_FI_gut, lambda6_FI_gut = mphi_run(phi, phigut, mphi_FI_phi0), lambda6_run(phi, phigut, lambda6_FI_phi0)
    if return_Acal:
        Acal_CORR = Acal * fac_corr
        return Acal_CORR
    if return_alpha:
        return Acal**2/A6**2
    return mphi_FI_gut, lambda6_FI_gut

def V(phi, phi0, A6gut, dv, infl_type):
    mphi_FI_gut, lambda6_FI_gut = FI(phi0, A6gut, dv, infl_type)
    return V_MSSM(phi, infl_type, mphi_FI_gut, A6gut, lambda6_FI_gut)


i = 33
infl_type = 'lle'
file_i = pd.read_csv('save_file2/lle_bp2_0_3047_9665.csv', dtype = str).iloc[i]
Vrge = lambda phi : V_MSSM(phi, infl_type, mp.mpf(file_i['mphi_gut']), mp.mpf(file_i['A6_gut']), mp.mpf(file_i['lambda6_gut']))
#print(mp.mpf(file_i['mphi_phi0']), mp.mpf(file_i['A6_phi0']), mp.mpf(file_i['lambda6_phi0']), mp.mpf(file_i['phi0']))
mp.plot([lambda phi : mp.log10(eps1_(Vrge, phi))], points = 200, xlim = (mp.mpf(file_i['phi0'])*0.9, mp.mpf(file_i['phi0'])*1.1))

"""
infl_type = 'lle'
i =33
file_i = pd.read_csv('lle_bp2_0_3047_9665.csv', dtype = str).iloc[i]
#Vrge = lambda phi :V(phi, mp.mpf(file_i['phi0']), mp.mpf(file_i['A6_phi0']), mp.mpf(file_i['dv']), infl_type)
#mp.plot([lambda phi : mp.log10(eps1_(Vrge, phi))], points = 200, xlim = (mp.mpf(file_i['phi0'])*0.9, mp.mpf(file_i['phi0'])*1.1))

#mphi_phi0, lambda6_phi0, phi0 = FI(mp.mpf(file_i['phi0']),
#                                   A6_run_lle_(mp.mpf(file_i['phi0']), mp.mpf('3e16'), mp.mpf(file_i['A6_phi0']),
#                                   m1gut, g1gut, m2gut, g2gut ), mp.mpf(file_i['dv']),infl_type,at_phi0=True)
#print(mphi_phi0, lambda6_phi0)
#print(mp.mpf(file_i['mphi_phi0']), mp.mpf(file_i['lambda6_phi0']))
Vrge = lambda phi :V_MSSM_run(phi, mp.mpf(file_i['phi0']),
                              infl_type, mp.mpf(file_i['mphi_phi0']), 
                              mp.mpf(file_i['A6_phi0']), mp.mpf(file_i['lambda6_phi0']))
#Vrge = lambda phi :V_MSSM_run(phi, mp.mpf(file_i['phi0']), infl_type, mphi_phi0, mp.mpf(file_i['A6_phi0']), lambda6_phi0)
mp.plot([lambda phi : mp.log10(eps1_(Vrge, phi))], points = 200, xlim = (mp.mpf(file_i['phi0'])*0.9, mp.mpf(file_i['phi0'])*1.1))
"""
"""
infl_type = 'lle'
A6 = mp.mpf('34890.76714837223375356279')
phi0 = mp.mpf('1000000000000000.0')
#A6 = mp.mpf('34614.18342539386515961698')
dv = mp.mpf('8306.090791034694650241906')
mphi, lambda6 = FI(phi0, A6, dv, infl_type)
print(mphi, A6, lambda6)
Vrge = lambda phi :V(phi, phi0, A6, dv, infl_type)
mp.plot([lambda phi : mp.log10(eps1_(Vrge, phi))], points = 200, xlim = (1e15*0.9, 1e15*1.1))
Vrge = lambda phi :V_MSSM(phi, infl_type, mphi, A6, lambda6)
mp.plot([lambda phi : mp.log10(eps1_(Vrge, phi))], points = 200, xlim = (1e15*0.9, 1e15*1.1))
"""

"""
A6 = mp.mpf('34890.76722504272163728756')
mphi = mp.mpf('7720.542680876346287167992')
lambda6 = mp.mpf('0.03189619730663271332928875')

A6 = mp.mpf('34614.18342539386515961698')
mphi = mp.mpf('7729.900012222098188527119')
lambda6 = mp.mpf('0.0351441983079566550744306')

Vrge = lambda phi :V_MSSM(phi, infl_type, mphi, A6, lambda6)
mp.plot([lambda phi : mp.log10(eps1_(Vrge, phi))], points = 200, xlim = (1e15*0.9, 1e15*1.1))
"""
"""
infl_type = 'lle'
mphi, A6, lambda6 = mp.mpf('7727.133865030157466888768'), mp.mpf('34773.79171405754952459827'), mp.mpf('0.03189368007223351785253577')
i = 20
file_i = pd.read_csv('save_file/lle_bp1_0_3047_9665.csv', dtype = str).iloc[i]
Vrge = lambda phi :V_MSSM(phi, infl_type, mp.mpf(file_i['mphi_gut']), mp.mpf(file_i['A6_gut']), mp.mpf(file_i['lambda6_gut']))
mp.plot([lambda phi : mp.log10(eps1_(Vrge, phi))], points = 200, xlim = (mp.mpf(file_i['phi0'])*0.9, mp.mpf(file_i['phi0'])*1.1))

i = 20
infl_type = 'tree'
file_i = pd.read_csv('tree_0_3047_9665.csv', dtype = str).iloc[i]
Vtree = lambda phi : V_MSSM(phi, infl_type, mp.mpf(file_i['mphi_phi0']), mp.mpf(file_i['A6_phi0']), mp.mpf(file_i['lambda6_phi0']))
print(mp.mpf(file_i['mphi_phi0']), mp.mpf(file_i['A6_phi0']), mp.mpf(file_i['lambda6_phi0']), mp.mpf(file_i['phi0']))
mp.plot([lambda phi : mp.log10(eps1_(Vtree, phi))], points = 200, xlim = (mp.mpf(file_i['phi0'])*0.9, mp.mpf(file_i['phi0'])*1.1))
"""
