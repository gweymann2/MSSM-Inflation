import mpmath as mp
from mpmath import *
import numpy as np
import pandas as pd
import time
import sys
import platform
import argparse

dps_def = 15
mp.dps = dps_def
extradps_def = 25
addprec_deriv = 30

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

def N_low(V, phimin, phimax):
    integ = mp.quad(lambda phi : V(phi)/DV(V, phi), [phimin, phimax], verbose=False, method='tanh-sinh')
#     print(nstr(-integ/Mp**2, 20), end=', ')
    return -integ/Mp**2

def find_phi_st(V, phi, phi_end, Pstar, lnRrad):
    kstar, lnMpcToKappa, HubbleSquareRootOf3OmegaRad, RelatDofRatio = mp.mpf('0.05'), mp.mpf('130.282'), mp.mpf('7.5437e-63'), mp.mpf('1')
    N0 = log(kstar) - lnMpcToKappa - 0.5*log(HubbleSquareRootOf3OmegaRad) -mp.mpf('0.25')*log(RelatDofRatio)
    Delta_N_star = N_low(V, phi, phi_end)
    return -Delta_N_star + lnRrad - N0 - 0.25*mp.log(9/eps1_(V, phi)/(3-eps1_(V, phi_end))*V(phi_end)/V(phi))+0.25*mp.log(8*mp.pi**2*Pstar)

def endinf(V, start):
    return mp.findroot(lambda phi: eps1_(V, phi)-1, 0.99*start)

def phi_star(V, Pstar, lnRrad, phi_single, plotose=False):
    phi_start = mp.findroot(lambda phi : D2V(V, phi), phi_single)
    phi_end = endinf(V, phi_start*0.98)
    kstar, lnMpcToKappa, HubbleSquareRootOf3OmegaRad, RelatDofRatio = mp.mpf('0.05'), mp.mpf('130.282'), mp.mpf('7.5437e-63'), mp.mpf('1')
    N0 = log(kstar) - lnMpcToKappa - 0.5*log(HubbleSquareRootOf3OmegaRad) -mp.mpf('0.25')*log(RelatDofRatio)
    if plotose == True:
        mp.plot([lambda phi : N_low(V, phi, phi_end), lambda phi : lnRrad - N0 - 0.25*mp.log(9/eps1_(V, phi)/(3-eps1_(V, phi_end))*V(phi_end)/V(phi))+0.25*mp.log(8*mp.pi**2*Pstar)], xlim = (phi_start*(1-shift_plot),phi_start*(1+shift_plot)), points=20)
        mp.plot([lambda phi : mp.log10(eps1_(V, phi))], xlim = (phi_start*(1-shift_plot),phi_start*(1+shift_plot)), points=100)
        mp.plot([lambda phi : eps2_(V, phi)], xlim = (phi_start*(1-shift_plot),phi_start*(1+shift_plot)), points=100)
    start_time = time.process_time()
    phi_star = mp.findroot(lambda phi : find_phi_st(V, phi, phi_end, Pstar, lnRrad), x0 = (phi_start-10,phi_start+10), verbose = False, method='ridder')
    print(int(time.process_time()-start_time),end='s ')
    return phi_star

def ns_(V, phi):
    return mp.fsub(mp.fsub(1,2*eps1_(V, phi)),eps2_(V, phi))

def P_star(V, phi_sta):
    Vstar = V(phi_sta)
    eps1star = eps1_(V, phi_sta)
    return Vstar/(Mp**4*24*mp.pi**2*eps1star)


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



def norm_potential(x,alpha,phi0):
    return x**mp.mpf('2')-mp.mpf('2')/3*alpha*x**6+alpha/mp.mpf('5')*x**10

def norm_eps1(x,alpha,phi0):
    return mp.mpf('450')*(mp.mpf('1')+alpha*x**4*(-mp.mpf('2')+x**4))**2/((phi0*x)**mp.mpf('2')*(mp.mpf('15')+alpha*x**4*(-mp.mpf('10')+mp.mpf('3')*x**4))**2)

def norm_eps2(x,alpha,phi0):
    return (mp.mpf('60')*(mp.mpf('15')+alpha*x**mp.mpf('4')*(mp.mpf('40')+x**mp.mpf('4')*(-mp.mpf('78')+alpha*(mp.mpf('20')+mp.mpf('3')*x**8)))))/(phi0**2*x**2*(mp.mpf('15')+alpha*x**mp.mpf('4')*(-mp.mpf('10')+mp.mpf('3')*x**4))**2)
       
def norm_eps3(x,alpha,phi0):
    return (mp.mpf('60')*(mp.mpf('1')+alpha*x**4*(-mp.mpf('2')+x**mp.mpf('4')))*(mp.mpf('225')+alpha*x**mp.mpf('4')*(-mp.mpf('1350')+x**4*(mp.mpf('3915')+alpha*(-mp.mpf('2100')+mp.mpf('20')*(mp.mpf('81')-mp.mpf('10')*alpha)*x**4+mp.mpf('15')*(-mp.mpf('99')+mp.mpf('20')*alpha)*x**mp.mpf('8')+mp.mpf('90')*alpha*x**mp.mpf('12')+mp.mpf('9')*alpha*x**16)))))/(x**2*(mp.mpf('15')+alpha*x**mp.mpf('4')*(-mp.mpf('10')+mp.mpf('3')*x**4))**2*(mp.mpf('15')+alpha*x**4*(mp.mpf('40')+x**4*(-mp.mpf('78')+alpha*(mp.mpf('20')+mp.mpf('3')*x**8)))))/phi0**2

def x_endinf(alpha,phi0):
    xstart = mp.mpf('0.9841521')
    return mp.findroot(lambda x : norm_eps1(x,alpha,phi0)-mp.mpf('1'), xstart, solver='halley', verbose = False)

def efold_primitive(x,alpha,phi0):
    aplus=-alpha+mp.sqrt((alpha**2-alpha)*mp.mpc('1','0'))
    aminus=-alpha-mp.sqrt((alpha**2-alpha)*mp.mpc('1','0'))
    bplus=mp.mpf('2')*(aplus+alpha/mp.mpf('3'))/(aplus-aminus)
    bminus=mp.mpf('2')*(aminus+alpha/mp.mpf('3'))/(aminus-aplus)
    return phi0**2*(mp.re(x**2/20+bplus/(10*mp.sqrt(aplus))*mp.atan(mp.sqrt(aplus)*x**2)+bminus/(10*mp.sqrt(aminus))*mp.atan(mp.sqrt(aminus)*x**2)))

def get_calfconst_rrad(lnRrad,Pstar,epsEnd,potEnd):
    cmbMeasuredQuart = mp.mpf('0.25')*mp.log(Pstar*mp.mpf('8')*mp.pi**2)
    kstar, lnMpcToKappa, HubbleSquareRootOf3OmegaRad, RelatDofRatio = mp.mpf('0.05'), mp.mpf('130.282'), mp.mpf('7.5437e-63'), mp.mpf('1')
    N0 = mp.log(kstar) - lnMpcToKappa - mp.mpf('0.5')*mp.log(HubbleSquareRootOf3OmegaRad) - mp.mpf('0.25')*mp.log(RelatDofRatio)
    return - N0 + cmbMeasuredQuart-mp.mpf('0.25')*mp.log(potEnd/(mp.mpf('3')-epsEnd)) + lnRrad

def x_rrad(alpha,phi0,lnRrad,Pstar):
    xEnd = x_endinf(alpha,phi0)
    epsOneEnd = norm_eps1(xEnd,alpha,phi0)
    potEnd = norm_potential(xEnd,alpha,phi0)
    primEnd = efold_primitive(xEnd,alpha,phi0)
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)
    calFplusNuEnd = calF+primEnd
    x_eps10 = 1
    return mp.findroot(lambda x : find_x_rrad(x, alpha, phi0,calFplusNuEnd), mp.mpf('1'), solver='halley', verbose = False)

def find_reheat_rrad_leadorder(nuStar,calFplusNuEnd,epsOneStar,Vstar):
    return nuStar - calFplusNuEnd + mp.mpf('0.25')*mp.log(mp.mpf('9')/(epsOneStar*Vstar))

def find_x_rrad(x,alpha,phi0,calFplusNuEnd):
    nuStar = efold_primitive(x,alpha,phi0)
    epsOneStar = norm_eps1(x,alpha,phi0)
    Vstar = norm_potential(x,alpha,phi0)
    res = find_reheat_rrad_leadorder(nuStar,calFplusNuEnd,epsOneStar,Vstar)
    return res

def ns_from_alpha(alpha, phi0B, lnRrad, Pstar):
    phi0 = phi0B*((5*alpha+mp.sqrt(mp.mpf('25')*alpha**2-mp.mpf('9')))/(9*alpha))**(-mp.mpf('0.25'))
    xstar = x_rrad(alpha, phi0, lnRrad, Pstar)
    return mp.mpf('1')-2*norm_eps1(xstar, alpha, phi0)-norm_eps2(xstar, alpha, phi0)

def alpha_coeff(phi0_B,lnRrad):
    phi0_B_Vec=[0.000001,0.000003,0.00001,0.0001,0.001,0.01,0.1]
    if lnRrad == 0:
      coeff_Vec=[2.5,2.2,2.05,1.7,1.4,1.2,1.]
    elif lnRrad == -10:
      coeff_Vec=[5.2,4.4,3.8,3.,2.4,2.,1.6]
    i=0
    while(phi0_B > phi0_B_Vec[i] and i < 5):
        i=i+1
    i = i-1
    return coeff_Vec[i]+(coeff_Vec[i+1]-coeff_Vec[i])*mp.log(phi0_B/phi0_B_Vec[i])/mp.log(phi0_B_Vec[i+1]/phi0_B_Vec[i])

def alpha_from_phi0B_and_ns(phi0B, ns, lnRrad,Pstar):
    alpha_min = mp.mpf('1')-alpha_coeff(phi0B,lnRrad)*phi0B**4*mp.pi**2/(mp.mpf('900')*mp.mpf('50')**2)
    alpha_max = mp.fsub(mp.mpf('1'),mp.mpf('1e-30'))
    return mp.findroot(lambda alpha : ns_from_alpha(alpha, phi0B, lnRrad, Pstar)-ns, (alpha_min, alpha_max), solver='anderson', verbose=False)

def aspic(lnRrad, ns_f, phi0_input, Pstar, _str = False):
    with extradps(30):
        phi0_input = [mp.re(phi0) for phi0 in phi0_input]
        phi0B_list, phistar_list, alpha_list, mphi_list, A6_list, lambda6_list, ns_list, r_list, alphas_list, eps10_list = [], [], [], [], [], [], [], [], [], []
        nphi0 = len(phi0_input)
        for i, phi0B in enumerate(phi0_input):

            alpha = alpha_from_phi0B_and_ns(phi0B/Mp, ns_f, lnRrad, Pstar)
            phi0 = phi0B*((5*alpha+mp.sqrt(mp.mpf('25')*alpha**2-mp.mpf('9')))/(9*alpha))**(-mp.mpf('0.25'))
            eps10 = norm_eps1(mp.mpf('1'),alpha,phi0/Mp)
            xstar = x_rrad(alpha, phi0/Mp, lnRrad, Pstar)
            eps1 = norm_eps1(xstar,alpha,phi0/Mp)
            eps2 = norm_eps2(xstar,alpha,phi0/Mp)
            eps3 = norm_eps3(xstar,alpha,phi0/Mp)

            ns, r, alpha_s = 1 - 2*eps1 - eps2, 16*eps1, -2*eps1*eps2-eps2*eps3
            M = (Pstar * 8 * mp.pi**2 * Mp**2 * eps1 * 3 * Mp**2 / norm_potential(xstar,alpha,phi0/Mp))**mp.mpf('0.25')
            mphiBoehm2 = 2 * M**4 / (phi0)**2 
            AB=mp.sqrt(80*alpha)*M**2/(phi0)
            lambdaB=Mp**3*mp.sqrt(alpha/5)*M**2/(phi0)**5

            if not _str:
                phi0B_list.append(phi0B)
                phistar_list.append(xstar*phi0)
                alpha_list.append(alpha)
                mphi_list.append(mphiBoehm2**0.5)
                A6_list.append(AB / fac_corr)            # !!!!!!!!! CORR !!!!!!!!! #
                lambda6_list.append(lambdaB)
                ns_list.append(ns)
                r_list.append(r)
                alphas_list.append(alpha_s)
                eps10_list.append(eps10)
            else:
                phi0B_list.append(nstr(phi0B, 32))
                phistar_list.append(nstr(xstar*phi0,32))
                alpha_list.append(nstr(alpha, 32))
                mphi_list.append(nstr(mphiBoehm2**0.5,32))
                A6_list.append(nstr(AB / fac_corr,32))            # !!!!!!!!! CORR !!!!!!!!! #
                lambda6_list.append(nstr(lambdaB,32))
                ns_list.append(nstr(ns,32))
                r_list.append(nstr(r,32))
                alphas_list.append(nstr(alpha_s, 32))
                eps10_list.append(nstr(eps10,32))
            
        return phi0B_list, phistar_list, alpha_list, mphi_list, A6_list, lambda6_list, ns_list, r_list, alphas_list, eps10_list


def find_dv(A6gut, phi0, dv, key, infl_type):
    mphi_FI_gut, lambda6_FI_gut = FI(phi0, A6gut, dv, infl_type)
    V_ = lambda phi : V_MSSM(phi, infl_type, mphi_FI_gut, A6gut, lambda6_FI_gut)
    phi_sta = phi_star(V_, mp.mpf('2.0989031673191437e-9'), mp.mpf('0'), phi0)
    global ns_star, P_sta, phi_sta_stock
    phi_sta_stock = phi_sta
    ns_star = ns_(V_, phi_sta)
    print(nstr(ns_star, 7), end = ', ')
    P_sta = P_star(V_, phi_sta)
    return ns_star

def dv_ns(A6gut, phi0, ns, infl_type, key='A6'):
    dv_start = mp.sqrt(2*aspic(mp.mpf('0'), ns, [phi0], mp.mpf('2.0989031673191437e-9'))[-1][0]*V(phi0, phi0, A6gut, mp.mpf('0'), infl_type)**2/Mp**2)
    t = time.process_time()
    dv_ns_res = mp.findroot(lambda dv : find_dv(A6gut, phi0, dv, key, infl_type) - ns, method='ridder', x0 = (dv_start*0.99, dv_start*1.01), tol=1e-6, verbose=False)
    global A6_stock, phi0_stock, dv_stock
    A6_stock, phi0_stock, dv_stock = A6gut, phi0, dv_ns_res
    print('     ', int(time.process_time()-t),'s', nstr(P_sta, 5))
    return P_sta

def A6_start(phi0, As, ns, infl_type):
    t=time.process_time()
    A6_start_start = mp.mpf('10000') # should be robust                      # udd       ! care here !
    Vphi0 = lambda A6gut : V(phi0, phi0, A6gut, mp.mpf('0'), infl_type)
    eps1starphi0 = lambda A6gut : aspic(mp.mpf('0'), ns, [phi0], As)[7][0]/16
    #mp.plot(lambda A6gut : Vphi0(A6gut)/(24*mp.pi**2*eps1starphi0(A6gut)*As*Mp**4)-1, xlim = (0,1000), points = 20)    
    try : A6_start_ = mp.findroot(lambda A6gut : Vphi0(A6gut)/(24*mp.pi**2*eps1starphi0(A6gut)*As*Mp**4)-1, x0 = (A6_start_start*0.9, A6_start_start*1.1), method='ridder', verbose=False)
    except : A6_start_ = mp.mpf('400')
    print(str(int(time.process_time()-t))+'s start,')
    return A6_start_

def A6_As(phi0, As, ns, infl_type):
    A6_start_ = A6_start(phi0, As, ns, infl_type)
    A6_As_res = mp.findroot(lambda A6 : dv_ns(A6, phi0, ns, infl_type, 'A6') - As, method='ridder', x0 = (A6_start_*0.99, A6_start_*1.01), tol=1e-6, verbose=False)
    return phi0_stock, A6_stock, dv_stock, phi_sta_stock

    
def phi0_start(A6gut, As, ns, infl_type):
    t=time.process_time()
    Vphi0 = lambda phi0 : V(phi0, phi0, A6gut, mp.mpf('0'), infl_type)
    eps1starphi0 = lambda phi0 : aspic(mp.mpf('0'), ns, [phi0], As)[7][0]/16
    print(str(int(time.process_time()-t))+'s start,')
    try: 
        phi0_start_start = mp.mpf('1e15')
        phi0_start = mp.findroot(lambda phi0 : (Vphi0(phi0)-24*mp.pi**2*eps1starphi0(phi0)*As*Mp**4)/Vphi0(phi0), (phi0_start_start*0.999, phi0_start_start*1.001), method='ridder', tol=1e-16, verbose=False)
    except:
        phi0_start_start = mp.mpf('1e16')
        phi0_start = mp.findroot(lambda phi0 : (Vphi0(phi0)-24*mp.pi**2*eps1starphi0(phi0)*As*Mp**4)/Vphi0(phi0), (phi0_start_start*0.999, phi0_start_start*1.001), method='ridder', tol=1e-16, verbose=False)
    return phi0_start


def phi0_As(A6gut, As, ns, infl_type):
    phi0_start_ = phi0_start(A6gut, As, ns, infl_type)
    print('ok', end='')
    phi0_As_res = mp.findroot(lambda phi0 : dv_ns(A6gut, phi0, ns, infl_type, 'phi0') - As, method='ridder', x0 = (phi0_start_*0.99, phi0_start_*1.01), tol=1e-6, verbose=False)
    return phi0_stock, A6_stock, dv_stock, phi_sta_stock

def contour_generator(phi0_list, A6_list, As, ns, lnRrad, infl_type, name_file): #### !!!! lnRrad and As only work at tree for now !!!! #####
    if A6_list.size == 0:
        print('Horizontal: ns search playing on dv')
        print('Vertical: As search playing on A6\n\n--------------------------------------------------\n')
        A6_list, dv_list, mphi_list, lambda6_list, phi_sta_list, alpha_list = [], [], [], [], [], []
        if infl_type=='tree':
            phi0_list_, phi_sta_list, alpha_list, mphi_list, A6_list, lambda6_list, ns_list, r_list, alphas_list, eps10_list = aspic(lnRrad, ns, phi0_list, As, _str = True)
            for i, phi0 in enumerate(phi0_list_):
                #print(eps10_list[i],mp.sqrt(mp.mpf('2')),mphi_list[i], A6_list[i], lambda6_list[i])#,V_MSSM(phi0, 'tree', mphi_list[i], A6_list[i], lambda6_list[i])/Mp)
                dv_list.append(mp.mpf(eps10_list[i])*mp.sqrt(mp.mpf('2'))*V_MSSM(mp.mpf(phi0), 'tree', mp.mpf(mphi_list[i]), mp.mpf(A6_list[i]), mp.mpf(lambda6_list[i]))/Mp)
        else:
            for k, phi0 in enumerate(phi0_list):
                print('Point '+str(k+1)+'/'+str(len(phi0_list)))
                print('phi0 = ',nstr(phi0,25))
                print('As = ',nstr(As,10))
                print('ns = ',nstr(ns,10),'\n')
                start_time = time.process_time() 
                phi0, A6, dv, phi_sta = A6_As(phi0, As, ns, infl_type)
    #            phi0, A6, dv, phi_sta = phi0_As(A6, As, ns, infl_type)

                time_step = time.process_time() - start_time
                print('\nPoint done in '+str(int(time_step)) +" s.")
                with extradps(16):
                    mphi_phi0, lambda6_phi0, A6_phi0 = FI(phi0, A6, dv, infl_type, at_phi0=True)
                    mphi, lambda6 = FI(phi0, A6, dv, infl_type)
                    alpha = FI(phi0, A6, 0, infl_type, return_Acal = True)**2/(1*fac_corr**2*40 * FI(phi0, A6, dv, infl_type, at_phi0=True)[0]**2)
    #                alpha_2 = FI(phi0, A6, dv, infl_type, return_alpha = True)
                    A6_list.append(nstr(A6_phi0,32))
                    dv_list.append(nstr(dv,32))
                    mphi_list.append(nstr(mphi_phi0,32))
                    lambda6_list.append(nstr(lambda6_phi0,32))
                    phi_sta_list.append(nstr(phi_sta,32))
                    alpha_list.append(nstr(alpha,32))
                    print('phi0 = ' + nstr(phi0, 25))
                    print('A6 = ' + nstr(A6_phi0, 25))
                    print('dv = ' + nstr(dv, 25))
                    print('mphi = ' + nstr(mphi_phi0, 25))
                    print('lambda6 = ' + nstr(lambda6_phi0, 25))
                    print('A6gut = ' + nstr(A6, 25))
                    print('mphigut = ' + nstr(mphi, 25))
                    print('lambda6gut = ' + nstr(lambda6, 25))
                    print('phi* = ' + nstr(phi_sta, 25))
                    print('alpha = ' + nstr(alpha, 25))
    #                print('alpha2 = ' + nstr(alpha_2, 25))
                    print('r = '+nstr(16*eps1_(lambda phi : V_MSSM(phi, infl_type, mphi, A6, lambda6), phi_sta),25))
                    print('\n--------------------------------------------------\n')
        data = pd.DataFrame(data={'phi0':[nstr(phi0_,32) for phi0_ in phi0_list], 'mphi_phi0':mphi_list, 'A6_phi0':A6_list, 'lambda6_phi0':lambda6_list, 'dv':dv_list, 'phi*':phi_sta_list, 'alpha':alpha_list}, dtype=str)
        data.to_csv(name_file)
        return name_file+' written.'
    elif phi0_list.size == 0:
        print('Horizontal: ns search playing on dv')
        print('Vertical: As search playing on A6\n\n--------------------------------------------------\n')
        phi0_list, dv_list, mphi_list, lambda6_list, phi_sta_list, alpha_list = [], [], [], [], [], []

        for k, A6 in enumerate(A6_list):
            print('Point '+str(k+1)+'/'+str(len(A6_list)))
            print('A6 = ',nstr(A6,25))
            print('As = ',nstr(As,10))
            print('ns = ',nstr(ns,10),'\n')
            start_time = time.process_time() 
            phi0, A6, dv, phi_sta = phi0_As(A6, As, ns, infl_type)

            time_step = time.process_time() - start_time
            print('\nPoint done in '+str(int(time_step)) +" s.")
            with extradps(16):
                mphi_phi0, lambda6_phi0, A6_phi0 = FI(phi0, A6, dv, infl_type, at_phi0=True)
                mphi, lambda6 = FI(phi0, A6, dv, infl_type)
                alpha = FI(phi0, A6, 0, infl_type, return_Acal = True)**2/(1*fac_corr**2*40 * FI(phi0, A6, dv, infl_type, at_phi0=True)[0]**2)
    #                alpha_2 = FI(phi0, A6, dv, infl_type, return_alpha = True)
                phi0_list.append(nstr(phi0,32))
                dv_list.append(nstr(dv,32))
                mphi_list.append(nstr(mphi_phi0,32))
                lambda6_list.append(nstr(lambda6_phi0,32))
                phi_sta_list.append(nstr(phi_sta,32))
                alpha_list.append(nstr(alpha,32))
                print('phi0 = ' + nstr(phi0, 25))
                print('A6 = ' + nstr(A6_phi0, 25))
                print('dv = ' + nstr(dv, 25))
                print('mphi = ' + nstr(mphi_phi0, 25))
                print('lambda6 = ' + nstr(lambda6_phi0, 25))
                print('A6gut = ' + nstr(A6, 25))
                print('mphigut = ' + nstr(mphi, 25))
                print('lambda6gut = ' + nstr(lambda6, 25))
                print('phi* = ' + nstr(phi_sta, 25))
                print('alpha = ' + nstr(alpha, 25))
    #                print('alpha2 = ' + nstr(alpha_2, 25))
                print('r = '+nstr(16*eps1_(lambda phi : V_MSSM(phi, infl_type, mphi, A6, lambda6), phi_sta),25))
                print('\n--------------------------------------------------\n')
        data = pd.DataFrame(data={'phi0':phi0_list, 'mphi_phi0':mphi_list, 'A6_phi0':A6_list, 'lambda6_phi0':lambda6_list, 'dv':dv_list, 'phi*':phi_sta_list, 'alpha':alpha_list}, dtype=str)
        data.to_csv(name_file)
        return name_file+' written.'


def name_file(infl_type, bp, lnRrad, ln1010As, ns, add=''):
    return ((infl_type=='tree')*infl_type)+((infl_type!='tree')*(infl_type+'_bp'+str(bp)))+'_'+str(int(-lnRrad))+'_'+str(int(ln1010As*1000))+'_'+str(int(10000*ns))+add+'.csv'

if __name__ == "__main__":

    print('\nall packages are successfully loaded')
    print('python version: ' + sys.version)
    print('Current OS: ' + platform.platform()+'\n')

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

        
    parser = argparse.ArgumentParser()
    parser.add_argument("-ns", "--ns", help= "ns", type = float)
    parser.add_argument("-ln1010As", "--ln1010As", help= "ln1010As", type = float)
    parser.add_argument("-lnRrad", "--lnRrad", help= "lnRrad", type = int)
    parser.add_argument("-infl_type", "--infl_type", help= "infl_type", type = str)
    parser.add_argument("-bp", "--bp", help= "bp", type = int)
    parser.add_argument("-xaxis", "--xaxis", help= "xaxis", type = str)
    parser.add_argument("-logorlin", "--logorlin", help= "logorlin", type = str)
    args = parser.parse_args()
    
    if args.ns == None: ns = mp.mpf('0.9665')
    else: ns = mp.mpf(str(args.ns))

    if args.ln1010As == None: ln1010As = 3.047
    else: ln1010As = args.ln1010As 
    As = mp.exp(mp.mpf(str(ln1010As)))*mp.mpf('10')**(-10)
    
    if args.lnRrad == None: lnRrad = mp.mpf('0')
    else: lnRrad = mp.mpf(str(args.lnRrad))
    
    if args.infl_type == None: infl_type = 'tree'
    else: infl_type = args.infl_type
    

    if args.bp == None and infl_type != 'tree': bp = 1
    else: bp = args.bp

    
    if args.xaxis == None: xaxis = 'phi0'
    else: xaxis = args.xaxis ; print('oj')

    if args.logorlin == None: logorlin = 'log'
    else: logorlin = args.logorlin
    
    if bp == 2 or bp == None:
        with extradps(extradps_def):
            g1gut = mp.mpf('0.704555660557172')
            g2gut = mp.mpf('0.690364970285155')
            g3gut = mp.mpf('0.684720653567032')
            m1gut = mp.mpf('897.774785812765')
            m2gut = mp.mpf('1789.66021839594')
            m3gut = mp.mpf('882.522487633969')
    elif bp == 1: 
        with extradps(extradps_def):
            g1gut = mp.sqrt(mp.mpf('5')/mp.mpf('3'))*mp.mpf('5.45185741e-01')
            g2gut = mp.mpf('6.90473022e-01')
            g3gut = mp.mpf('6.84972506e-01')
            m1gut = mp.mpf('1.36108022e+02')
            m2gut = mp.mpf('1.14286222e+03')
            m3gut = mp.mpf('8.98639714e+02')

    
    phi0_list0 = [mp.mpf(str(x)) for x in np.logspace(np.log10(100000000000000), np.log10(170000000000000), 16)][:-1]
    phi0_list1 = [mp.mpf(str(x)) for x in np.logspace(np.log10(170000000000000), np.log10(500000000000000), 31)][:-1]
    phi0_list2 = [mp.mpf(str(x)) for x in np.logspace(np.log10(500000000000000), np.log10(3e16), 15)]
    phi0_list_all = np.array(phi0_list0+phi0_list1+phi0_list2)
    phi0_list_rge = np.array(phi0_list1+phi0_list2)
    phi0_list_lin = np.linspace(1e14, 3e16, 100)

    #A6_list = np.array([40000])
    #A6_list = np.logspace(np.log10(1000), np.log10(2e8), 30)    
    A6_list = np.logspace(np.log10(3000), np.log10(2e8), 30)
    #A6_list = np.array([3000])
    #A6_list = np.logspace(np.log10(6000), np.log10(2e7), 30)

    add = ''
    if infl_type == 'tree': phi0_list = phi0_list_all
    else : phi0_list = phi0_list_rge
    if logorlin == 'lin': 
        phi0_list = phi0_list_lin
        add += '_lin'
        
    if xaxis == 'A6': 
        phi0_list, A6_list = np.array([]), A6_list
        add += '_A6'
    elif xaxis == 'phi0': 
        phi0_list, A6_list  = phi0_list, np.array([])

    # remove
    #A6_list = np.array([]) 
    #phi0_list = np.array([1.2e15]) 
    #add = '_phi012'
    #
    phi0_list = [mp.mpf('1202271120486381.0')]#[mp.mpf('1e15')]
    add = '_single'
    _name_file = name_file(infl_type, bp, lnRrad, ln1010As, ns, add = add)
    #
    #_name_file = 'xaxis_A6_long/'+_name_file
    #
    print('job launched with ns = '+ str(ns)+ ', '+ 'As = '+ str(ln1010As)+ ', '+ 'lnRrad = '+ str(lnRrad)+ ', '+ 'infl_type = '+ infl_type+', add ='+add+' =====> '+_name_file+'\n')    
    #print(phi0_list,A6_list,)

    #print(m1gut, m2gut, )
    #phi0 = mp.mpf('1000000000000000.0')
    #A6 = mp.mpf('34890.76714837223375356279')
    #dv = mp.mpf('8306.090791432641501792549') 
    #print(FI(phi0, A6, dv, 'lle'))
    contour_generator(phi0_list,
                      A6_list,
                      As,
                      ns,
                      lnRrad,
                      infl_type,
                      _name_file,
                      ) 

