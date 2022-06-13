import mpmath as mp
import os
from mpmath import *
import numpy as np
#import pandas as pd
import time
import matplotlib.pyplot as plt
import suspect_bind
import tree_level
import scipy.optimize as optimize

# put here a file at any scale (eg EWSB but GUT should work aswell) 
#initializing_file = "InflationNew_udd_hFun.in"
initializing_file = "InflationNew_udd_CoAnn.in"
#initializing_file = "InflationNew_udd_AFun.in"

aa=suspect_bind.suspect()
m_SLHA=aa.SLHAblock()
aa.Initialize(initializing_file)
aa.Execute()
#aa.addScale(2713)
#print(m_SLHA.EXTPAR(44, 2713))


def addscale_and_run(phi0):
    aa.addScale(phi0)
#    global g10, g20, g30, M10, M20, M30, A60
    g10, g20, g30 = m_SLHA.GAUGE(1,phi0)*np.sqrt(5/3),m_SLHA.GAUGE(2,phi0),m_SLHA.GAUGE(3,phi0)
    M10,M20,M30=m_SLHA.MSOFT(1,phi0),m_SLHA.MSOFT(2,phi0),m_SLHA.MSOFT(3,phi0)
    A60 = m_SLHA.INFLATION(1,phi0)
    return g10, g20, g30, M10, M20, M30, A60


phigut = 3e16
mp.dps = 16
lnMpinGev = 42.334
Mp = np.exp(lnMpinGev)
pre = 1/(8*np.pi**2)
b1, b2, b3 = 33/5, 1, -3
"""
g1gut = np.sqrt(5/3)*5.45185741e-01
g2gut = 6.90473022e-01
g3gut = 6.84972506e-01
m1gut = 1.36108022e+02
m2gut = 1.14286222e+03
m3gut = 8.98639714e+02
A6gut = 40000

g1 = lambda phi : g1gut/(np.sqrt(1-pre*b1*g1gut**2*mp.log(phi/phigut)))
g2 = lambda phi : g2gut/(np.sqrt(1-pre*b2*g2gut**2*mp.log(phi/phigut)))
g3 = lambda phi : g3gut/(np.sqrt(1-pre*b3*g3gut**2*mp.log(phi/phigut)))
M1 = lambda phi : m1gut*(g1(phi)/g1gut)**2
M2 = lambda phi : m2gut*(g2(phi)/g2gut)**2
M3 = lambda phi : m3gut*(g3(phi)/g3gut)**2

mphi_run_lle = lambda phi_start, phi_end, mphi_start : np.sqrt(mphi_start**2+(M2(phi_start)**2-M2(phi_end)**2)+1/11*(M1(phi_start)**2-M1(phi_end)**2))
A6_run_lle = lambda phi_start, phi_end, A6_start : A6_start-6*(M2(phi_start)-M2(phi_end))-6/11*(M1(phi_start)-M1(phi_end))
lambda6_run_lle = lambda phi_start, phi_end, lambda6_start : lambda6_start*(g2(phi_start)/g2(phi_end))**6*(g1(phi_start)/g1(phi_end))**(6/11)
mphi_run_udd = lambda phi_start, phi_end, mphi_start : mp.sqrt(mphi_start**2-mp.mpf('8')/9*(M3(phi_start)**2-M3(phi_end)**2)+mp.mpf('4')/99*(M1(phi_start)**2-M1(phi_end)**2))
A6_run_udd = lambda phi_start, phi_end, A6_start : A6_start+16/3*(M3(phi_start)-M3(phi_end))-8/33*(M1(phi_start)-M1(phi_end))
lambda6_run_udd = lambda phi_start, phi_end, lambda6_start : lambda6_start*(g3(phi_start)/g3(phi_end))**(mp.mpf('6')*(-mp.mpf('8'))/9)*(g1(phi_start)/g1(phi_end))**(4*mp.mpf('6')/99)
mphi_run_tree = lambda phi_start, phi_end, mphi_start : mphi_start
A6_run_tree = lambda phi_start, phi_end, mphi_start : A6_start
lambda6_run_tree = lambda phi_start, phi_end, lambda6_start : lambda6_start
"""

def epsilon_from_phi(phi, phi0):
    return phi/phi0-1

def phi_from_epsilon(epsilon, phi0):
    return phi0*(1+epsilon)

    
def V_and_Vprime(phi0, nu, epsilon, infl_type, g10,g20,g30,M10,M20,M30,A60, printt=False, give_param=False, suspect=True):
  
    if infl_type == 'lle':
        ga, gb, Ma, Mb, A6 = g10, g20, M10, M20, A60
        Ba, Bb, ea, eb = b1/(16*np.pi**2), b2/(16*np.pi**2), 1/11*b1/(4*np.pi**2), b2/(4*np.pi**2)

    elif infl_type == 'udd':
        ga, gb, Ma, Mb, A6 = g10, g30, M10, M30, A60
        Ba, Bb, ea, eb = b1/(16*np.pi**2), b3/(16*np.pi**2), 1/(15*np.pi**2), 2/(3*np.pi**2),# 4/99*b1/(4*np.pi**2), -8/9*b3/(4*np.pi**2)

    elif infl_type == 'tree':
        ga, gb, Ma, Mb, A6 = 0, 0, 0, 0, A60
        Ba, Bb, ea, eb = 0, 0, 0, 0
        
#        elif infl_type == 'yukawas':
#            ga, gb, gc, Ma, Mb, Mc, A6, Yukawas ... = g10, g20, g30, M10, M20, M30, A60, Yukawas0 ...            

    else: print('this infl_type is not implemented')

#    print(ga, gb, Ma, Mb, A6)

    global nbr_call
    nbr_call += 1
    
    phi = phi_from_epsilon(epsilon, phi0) 

    Ca=2*Ba
    Cb=2*Bb
    Fa=3*ea
    Fb=3*eb
    Da=Fa/2
    Db=Fb/2

    betaga = 2*Ba*ga**4
    betagb = 2*Bb*gb**4
    betaMa = Ca*ga**2*Ma
    betaMb = Cb*gb**2*Mb
    betam = -ea*ga**2*Ma**2-eb*gb**2*Mb**2
    betaA = Fa*ga**2*Ma+Fb*gb**2*Mb
    betalambda = -(Da*ga**2+Db*gb**2)
    
    CST1p = 2*nu/phi0-betam
    CST2p = -2*((Ba+Ca)*ea*Ma**2*ga**4+(Bb+Cb)*eb*Mb**2*gb**4)+3*betam 
#     print(CST1p, CST2p)
    c1_ = 12*(5+betalambda) 
    c2_ = 12*(Da*betaga+Db*betagb-(5+betalambda)*(9+2*betalambda)) 
    c3_ = betaA+A6*(6+betalambda) 
    c4_ = -(2*Ba+Ca)*Fa*Ma*ga**4-(2*Bb+Cb)*Fb*Mb*gb**4-betaA*(11+2*betalambda)+A6*(Da*betaga+Db*betagb-(5+betalambda)*(6+betalambda)) 
    
    COEF1 = 3*(CST2p*c1_ - CST1p*c2_)**2 + (CST2p*c3_ - CST1p*c4_)*(-c2_*c3_ + c1_*c4_) 
    COEF2 = 1/10*(6*(c1_ + c2_)*(CST2p*c1_ - CST1p*c2_) + (c3_ + c4_)*(-c2_*c3_ + c1_*c4_))
    COEF3 = 3/100*(c1_ + c2_)**2 
    r1 = COEF1/COEF3
    r2 = COEF2/COEF3
    Acal = np.sqrt(-r2+np.sqrt(-4*r1+r2**2))
    mphi = Acal/np.sqrt(40)
    lambda6 = 3*Mp**3/phi0**4*(CST2p*c1_ - CST1p*c2_ + 2*mphi**2*(c1_ + c2_))/(c2_*c3_ - c1_*c4_)
    V0 = 0.5*mphi**2*phi0**2-lambda6*A6/(6*Mp**3)*phi0**6+lambda6**2*phi0**10/Mp**6
    
    if give_param:
        return mphi, lambda6
    
#     mphiphi_func, A6phi_func, lambda6phi_func = lambda phi : float(mphi_run(phi0, phi, mphi)), lambda phi : float(A6_run(phi0, phi, A6)), lambda phi : float(lambda6_run(phi0, phi, lambda6))
    
#     V_func = lambda phi : 0.5*mphiphi_func(phi)**2*phi**2-lambda6phi_func(phi)*A6phi_func(phi)/(6*Mp**3)*phi**6+lambda6phi_func(phi)**2*phi**10/Mp**6
#     V = V_func(phi)
    
    x0 = phi0**4*lambda6/Mp**3
    c0 = nu/phi0
    c1 = -c0
    c2 = ea*(-309*x0**2+0.5*Ma**2+27/4*(A6-2*Ma)*x0)*ga**2+eb*(-309*x0**2+0.5*Mb**2+27/4*(A6-2*Mb)*x0)*gb**2-6*A6*x0+280*x0**2
    c3 = ea*(-898*x0**2-1/3*Ma**2+15/2*(A6-2*Ma)*x0)*ga**2+eb*(-898*x0**2-1/3*Mb**2+15/2*(A6-2*Mb)*x0)*gb**2-4*A6*x0+560*x0**2
    c4 = ea*(-3085/2*x0**2+1/4*Ma**2+27/8*(A6-2*Ma)*x0)*ga**2+eb*(-3085/2*x0**2+1/4*Mb**2+27/8*(A6-2*Mb)*x0)*gb**2-A6*x0+700*x0**2
    c5 = ea*(-1654*x0**2-1/5*Ma**2+3/10*(A6-2*Ma)*x0)*ga**2+eb*(-1654*x0**2-1/5*Mb**2+3/10*(A6-2*Mb)*x0)*gb**2+560*x0**2
    c6 = ea*(-1107*x0**2+1/6*Ma**2-1/20*(A6-2*Ma)*x0)*ga**2+eb*(-1107*x0**2+1/6*Mb**2-1/20*(A6-2*Mb)*x0)*gb**2+280*x0**2
    c7 = ea*(-3054/7*x0**2-1/7*Ma**2+1/70*(A6-2*Ma)*x0)*ga**2+eb*(-3054/7*x0**2-1/7*Mb**2+1/70*(A6-2*Mb)*x0)*gb**2+80*x0**2
    c8 = ea*(-2367/28*x0**2+1/8*Ma**2-3/560*(A6-2*Ma)*x0)*ga**2+eb*(-2367/28*x0**2+1/8*Mb**2-3/560*(A6-2*Mb)*x0)*gb**2+80*x0**2
    c9 = ea*(-10/3*x0**2-1/9*Ma**2+1/420*(A6-2*Ma)*x0)*ga**2+eb*(-10/3*x0**2-1/9*Mb**2+1/420*(A6-2*Mb)*x0)*gb**2
    c10 = ea*(1/3*x0**2+1/10*Ma**2-1/840*(A6-2*Ma)*x0)*ga**2+eb*(1/3*x0**2+1/10*Mb**2-1/840*(A6-2*Mb)*x0)*gb**2
    
    Vbis = V0+phi0**2*(epsilon * c0 + epsilon**2/2*(c0+c1) + epsilon**3/3*(c1+c2) + epsilon**4/4*(c2+c3) + epsilon**5/5*(c3+c4) + epsilon**6/6*(c4+c5) + epsilon**7/7*(c5+c6) + epsilon**8/8*(c6+c7) + epsilon**9/9*(c7+c8) + epsilon**10/10*(c8+c9) + epsilon**11/11*(c9+c10))
    Vprime = phi*(c0+c1*epsilon+epsilon**2*c2+epsilon**3*c3+epsilon**4*c4+epsilon**5*c5+epsilon**6*c6+epsilon**7*c7+epsilon**8*c8+epsilon**9*c9+epsilon**10*c10)
#     print(Vprime)
    Vprimeprime = 2*epsilon*(-nu/phi0+c2)+3*epsilon**2*(c2+c3)+4*epsilon**3*(c3+c4)+5*epsilon**4*(c4+c5)+6*epsilon**5*(c5+c6)+7*epsilon**6*(c6+c7)+8*epsilon**7*(c7+c8)+9*epsilon**8*(c8+c9)+10*epsilon**9*(c9+c10)+11*epsilon**10*c10
    
    epsilonend = epsilon_from_phi(900809386587443.9, phi0)
    a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 = V0, phi0**2*c0, phi0**2*(c0+c1)/2, phi0**2*(c1+c2)/3, phi0**2*(c2+c3)/4, phi0**2*(c3+c4)/5, phi0**2*(c4+c5)/6, phi0**2*(c5+c6)/7, phi0**2*(c6+c7)/8, phi0**2*(c7+c8)/9, phi0**2*(c8+c9)/10, phi0**2*(c9+c10)/11
#     print(1/phi0*(a1+ 2*a2*epsilon+ a3*3*epsilon**2+ a4*4*epsilon**3+ a5*5*epsilon**4+ a6*6*epsilon**5+ a7*7*epsilon**6+ a8*8*epsilon**7+ a9*9*epsilon**8+ a10*10*epsilon**9+ a11*11*epsilon**10))
#     primVoverVprime_end = determine_intVoverVprime(epsilonend, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11)
#     primVoverVprime_phi = determine_intVoverVprime(epsilon, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11)
    L = [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11] 
#     intVoverVprime_phi_end = primVoverVprime_end - primVoverVprime_phi
    
    return float(Vbis), float(Vprime), float(Vbis/Vprime), float(Vprimeprime)#, L#, float(-1/Mp**2*intVoverVprime_phi_end)

nbr_call = 0

"""
def for_plotose(phi0, nu, infl_type, epsilonend):

    s = lambda phi : V_and_Vprime(phi0, nu, epsilon_from_phi(phi, phi0), infl_type)

    Vend, Vprime_end, VoverVprime_end, Vprimeprime_end = V_and_Vprime(phi0, nu, epsilonend, infl_type)
    epsilonmin = epsilonend
    #print(V_and_Vprime(phi0, nu, 0, infl_type)[1])
    N_left = lambda epsilonmax : float(mp.quad(lambda epsilon :  phi0*V_and_Vprime(phi0, nu, epsilon, infl_type)[2]/Mp**2, [epsilonmin, epsilonmax], verbose=False, method='tanh-sinh'))
    lnRrad, Pstar = 0, 2.0989031673191437e-9
    kstar, lnMpcToKappa, HubbleSquareRootOf3OmegaRad, RelatDofRatio = 0.05, 130.282, 7.5437e-63, 1
    N0 = mp.log(kstar) - lnMpcToKappa - 0.5 * mp.log(HubbleSquareRootOf3OmegaRad) - 0.25 * mp.log(RelatDofRatio)
    eps1phi = lambda phi : Mp**2/s(phi)[2]**2/2
    eps1phiend = Mp**2/VoverVprime_end**2/2
    N_right = lambda phi : lnRrad - N0 - 0.25*mp.log(9/eps1phi(phi)/(3-1)*Vend/s(phi)[0])+0.25*mp.log(8*np.pi**2*Pstar)
    #print(N_right(895378491085844.3883823681), N_left(epsilon_from_phi(895378491085844.3883823681,phi0)))

    ext = 80000
    phi_list1, phi_list2 = np.linspace(-ext,ext,20), np.linspace(-ext,ext,100)
    eps1list, eps2list, Nrightlist, Nleftlist = [], [], [], []
    for phimoinsphi0 in phi_list2:
        eps1list.append((lambda phimoinsphi0 : mp.log10(Mp**2/s(phimoinsphi0+phi0)[2]**2/2))(phimoinsphi0))
        eps2list.append((lambda phimoinsphi0 : 4*Mp**2/s(phimoinsphi0+phi0)[2]**2/2-2*Mp**2*(s(phimoinsphi0+phi0)[3]/s(phimoinsphi0+phi0)[0]))(phimoinsphi0))
    for phimoinsphi0 in phi_list1:
        Nrightlist.append((lambda phimoinsphi0 : N_right(phimoinsphi0+phi0))(phimoinsphi0))
        Nleftlist.append((lambda phimoinsphi0 : N_left(epsilon_from_phi(phimoinsphi0+phi0, phi0)))(phimoinsphi0))
    f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6,10), sharex=True)
    ax1.plot(phi_list1, Nrightlist, label='Late time')
    ax1.plot(phi_list1, Nleftlist, label='SR integral')
    ax1.legend(fontsize=13)
    ax1.set_ylabel(r'$N$',fontsize=13)
    ax1.grid()
    ax2.plot(phi_list2, eps1list)
    ax2.set_ylabel(r'$\epsilon_1$',fontsize=13)
    ax2.grid()
    ax3.plot(phi_list2, eps2list)
    ax3.set_xlabel(r'$\phi-\phi_0$',fontsize=13)
    ax3.set_ylabel(r'$\epsilon_2$',fontsize=13)
    ax3.grid()
    f.show()
"""

def find_phistar(phi, phi0, nu, infl_type, epsilonend, Vend, g10,g20,g30,M10,M20,M30,A60):
    epsilon = epsilon_from_phi(phi, phi0)
    V, Vprime, VoverVprime, Vprimeprime = V_and_Vprime(phi0, nu, epsilon, infl_type, g10,g20,g30,M10,M20,M30,A60)
    epsilonmin = epsilonend 
    N_left = float(mp.quad(lambda epsilon :  phi0*V_and_Vprime(phi0, nu, epsilon, infl_type, g10,g20,g30,M10,M20,M30,A60)[2]/Mp**2, [epsilonmin, epsilon_from_phi(phi, phi0)], verbose=False, method='tanh-sinh'))
    lnRrad, Pstar = 0, 2.0989031673191437e-9
    kstar, lnMpcToKappa, HubbleSquareRootOf3OmegaRad, RelatDofRatio = 0.05, 130.282, 7.5437e-63, 1
    N0 = mp.log(kstar) - lnMpcToKappa - 0.5 * mp.log(HubbleSquareRootOf3OmegaRad) - 0.25 * mp.log(RelatDofRatio)
    N_right = lnRrad - N0 - 0.25*mp.log(9/(Mp**2/VoverVprime**2/2)/(3-1)*Vend/V)+0.25*mp.log(8*np.pi**2*Pstar)
    #print(end='.')
    return N_right - N_left

def determine_phiend(phi, phi0, nu, infl_type, g10,g20,g30,M10,M20,M30,A60, verbose=False):
    find_phiend = lambda phi : Mp**2/V_and_Vprime(phi0, nu, epsilon_from_phi(phi, phi0), infl_type, g10,g20,g30,M10,M20,M30,A60)[2]**2/2-1
#     return mp.findroot(find_phiend, phi0*0.98, verbose = verbose)
    return optimize.ridder(find_phiend, phi0*0.99, phi0)

def newnew_phistar(phi0, nu, infl_type,g10, g20, g30, M10, M20, M30, A60, verbose=False, plotose=False):
 
    t=time.time()
    global nbr_call
    nbr_call = 0
    s = lambda phi : V_and_Vprime(phi0, nu, epsilon_from_phi(phi, phi0), infl_type, g10,g20,g30,M10,M20,M30,A60)
    phiend = determine_phiend(phi, phi0, nu, infl_type, g10,g20,g30,M10,M20,M30,A60)
    #print(phiend)
    epsilonend = epsilon_from_phi(phiend, phi0)
    Vend, Vprime_end, VoverVprime_end, Vprimeprime_end = V_and_Vprime(phi0, nu, epsilonend, infl_type, g10,g20,g30,M10,M20,M30,A60)
#     print(nbr_call)

    if plotose:
        for_plotose(phi0, nu, infl_type, epsilonend)
    
    find_phistar_phi = lambda phi : find_phistar(phi, phi0, nu, infl_type, epsilonend, Vend, g10,g20,g30,M10,M20,M30,A60)
    #phistar = optimize.ridder(find_phistar_phi, phi0*0.9999999, phi0)
#     phistar = optimize.brentq(find_phistar_phi, phi0*0.9999999, phi0)
    phistar = mp.findroot(find_phistar_phi, x0 = (phi0-10,phi0+10), verbose = verbose, method='ridder', tol = 1e-16)
#     print(nbr_call)

    sstar = s(phistar)
    eps1star = Mp**2/sstar[2]**2/2
    eps2star = 4*eps1star-2*Mp**2*(sstar[3]/sstar[0])
    As = sstar[0]/(Mp**4*24*np.pi**2*eps1star)
    ns = 1-2*eps1star-eps2star
#     print(nbr_call)

    return phistar, As, ns  

#nu = 3130.037484266968450061235
#phi0 = 895378491103795.0331033468
#infl_type = 'udd'
#print(V_and_Vprime(phi0, nu, 0, infl_type, *addscale_and_run(phi0)))

#phistar, As, ns = newnew_phistar(phi0, nu, infl_type, verbose=True, plotose=True)
#print(phistar, As, ns)

#print('phistar =', phistar)
#print('As =', "%.5e"%As)
#print('ns =', round(ns,6))



def find_dv(phi0, nu, key, infl_type, g10,g20,g30,M10,M20,M30,A6):
    global phistar_, As_, ns_
    global nbr_call
    phistar_, As_, ns_ = newnew_phistar(phi0, nu, infl_type, g10,g20,g30,M10,M20,M30,A6)
    print(nstr(mp.mpf(str(ns_)), 7), nbr_call, end = ', ')
    return ns_

def dv_ns(phi0, ns, As, infl_type, key='A6'):
    dv_start = tree_level.tree(phi0, As, ns, lnRrad)['dv']
    t = time.process_time()

    dv_ns_res = mp.findroot(lambda dv : find_dv(phi0, dv, key, infl_type, *addscale_and_run(phi0)) - ns, method='ridder', x0 = (dv_start*0.99, dv_start*1.01), tol=1e-6, verbose=False)
    global phi0_stock, dv_stock
    phi0_stock, dv_stock = phi0, dv_ns_res
    print('     ', int(time.process_time()-t),'s', nstr(mp.mpf(str(As_)), 5))
    return As_

def phi0_start(As, ns, infl_type):
    t=time.process_time()
    phi0_start_start = 1e15 # should be robust (more worried than above)
    eps1starphi0 = lambda phi0 : tree_level.tree(phi0, As, ns, 0)['eps10']
    Vphi0 = lambda phi0 : V_and_Vprime(phi0, 1, 0, infl_type,*addscale_and_run(phi0))[0]
    start = mp.findroot(lambda phi0 : (Vphi0(phi0)-24*np.pi**2*eps1starphi0(phi0)*As*Mp**4)/Vphi0(phi0), (phi0_start_start*0.9, phi0_start_start*1.1), method='ridder', tol=1e-30, verbose=False)
    print(str(int(time.process_time()-t))+'s start,')
    return float(start)

def phi0_As(As, ns, infl_type):
    phi0_start_ = phi0_start(As, ns, infl_type)
    phi0_As_res = mp.findroot(lambda phi0 : dv_ns(phi0, ns, As, infl_type, 'phi0') - As, method='ridder', x0 = (phi0_start_*0.99, phi0_start_*1.01), tol=1e-6, verbose=False)
    return phi0_stock, dv_stock, phistar_

#### comeback to gut and write out file
L = []
def mphi0_vs_msq1gut(msq1gut, phi0):

    init_slha = initializing_file
    os.system("cp "+init_slha+' '+output_file)
    msq1gut_str = str(msq1gut)
    slha_in = output_file

    msq1gut_str = str(msq1gut)
    with open(slha_in, 'r') as file:
        data = file.read()
        msq1_to_replace = data[3818:3825]
        print(msq1_to_replace)
        data = data.replace(msq1_to_replace, msq1gut_str)
        print(data)
        
    with open(slha_in, 'w') as file:
        file.write(data)
    #aa2 = suspect_bind.suspect()
    L.append(suspect_bind.suspect())
    m_SLHA=L[-1].SLHAblock()
    L[-1].Initialize(slha_in)
    L[-1].Execute()
    L[-1].addScale(phi0)
    mphi0 = m_SLHA.INFLATION(3,phi0)

    lambdaphi0 = m_SLHA.INFLATION(2,phi0)
    L[-1].addScale(3e16)
    lambda6gut = m_SLHA.INFLATION(2,3e16)
    global facteur, mphigut
    mphigut = m_SLHA.INFLATION(3,3e16)
    facteur = lambda6gut/lambdaphi0
    return mphi0



def msq1gut_vs_mphi0(mphi0_tilde, phi0):
    return mp.findroot(lambda msq1gut : mphi0_vs_msq1gut(msq1gut, phi0)-mphi0_tilde, mphi0_tilde, tol = 1e-10, verbose=False)


mp.dps = 16
infl_type = 'udd'
lnRrad = 0
As = 2.0989031673191437e-9
ns = 0.9653

phi0_res, dv_res, phi_sta_res = phi0_As(As, ns, infl_type)

#phi0_res = 545613800018766.60905
#dv_res = 41.544272492646398023

mphiphi0, lambda6phi0 = V_and_Vprime(phi0_res, dv_res, 0, infl_type,*addscale_and_run(phi0_res), give_param=True)
#output_file = 'being_modified.in'
#msq1gut_tilde = msq1gut_vs_mphi0(mphiphi0, phi0_res)
#init_gut = addscale_and_run(3e16)


print('\n--------------------------------------------------------------------------------------------------------------------------------')

print('\nIn input file '+str(initializing_file)+':')
print('infl_type\tg1_gut\t\tg2_gut\t\tg3_gut\t\tM1_gut\t\tM2_gut\t\tM3_gut\t\tA6_gut\n'+infl_type+'\t\t'+nstr(mp.mpf(str(init_gut[0])), 8)+'\t'+nstr(mp.mpf(str(init_gut[1])), 8)+'\t'+nstr(mp.mpf(str(init_gut[2])), 8)+'\t'+nstr(mp.mpf(str(init_gut[3])), 8)+'\t'+nstr(mp.mpf(str(init_gut[4])), 8)+'\t'+nstr(mp.mpf(str(init_gut[5])), 8)+'\t'+nstr(mp.mpf(str(init_gut[6])), 8))
print('\nCosmo constraints:')
print('As\t\tns\t\tlnRrad\n'+nstr(mp.mpf(As), 8)+'\t'+nstr(mp.mpf(str(ns)), 8)+'\t\t'+str(lnRrad))

print('\nDof tuning:')
print('phi0_res\tdv_res\n' + nstr(mp.mpf(str(phi0_res)), 8)+'\t'+nstr(mp.mpf(str(dv_res)), 8))

print('\nImplies for params at phi0:')
print('mphi_0\t\tlambda6_0\n' + nstr(mp.mpf(str(mphiphi0)), 8)+'\t'+nstr(mp.mpf(str(lambda6phi0)), 8))

print('to be run at GUT, which can be done through just below (needs to be improved)') 
#print('\nIn output file '+str(output_file)+':')
#print('msq1EWSB\tmphigut\t\tlambda6gut\n' + nstr(mp.mpf(str(msq1gut_tilde)), 8)+'\t'+nstr(mp.mpf(str(mphigut)), 8)+'\t'+nstr(mp.mpf(str(lambda6phi0*facteur)), 8))


print('\n--------------------------------------------------------------------------------------------------------------------------------\n')
