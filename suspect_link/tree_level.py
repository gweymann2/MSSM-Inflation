import mpmath as mp
from mpmath import *
def tree(phi0, As, ns, lnRrad):
    
    lnMpinGev = mp.mpf('42.334')
    Mp = mp.exp(lnMpinGev)
    
    def norm_potential(x,alpha,phi0):
        return x**mp.mpf('2')-mp.mpf('2')/3*alpha*x**6+alpha/mp.mpf('5')*x**10
    
    def norm_deriv(x,alpha,phi0):
        return mp.mpf('2')*(x-mp.mpf('2')*alpha*x**5+alpha*x**9)

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

    def aspic(lnRrad, ns_f, phi0B, Pstar):

            alpha = alpha_from_phi0B_and_ns(phi0B/Mp, ns_f, lnRrad, Pstar)
            phi0 = phi0B*((5*alpha+mp.sqrt(mp.mpf('25')*alpha**2-mp.mpf('9')))/(9*alpha))**(-mp.mpf('0.25'))
            xstar = x_rrad(alpha, phi0/Mp, lnRrad, Pstar)
            eps1 = norm_eps1(xstar,alpha,phi0/Mp)
            eps2 = norm_eps2(xstar,alpha,phi0/Mp)
            eps3 = norm_eps3(xstar,alpha,phi0/Mp)

            ns, r, alpha_s = 1 - 2*eps1 - eps2, 16*eps1, -2*eps1*eps2-eps2*eps3
            M = (Pstar * 8 * mp.pi**2 * Mp**2 * eps1 * 3 * Mp**2 / norm_potential(xstar,alpha,phi0/Mp))**mp.mpf('0.25')

            V0 = M**4*norm_potential(phi0B/phi0,alpha,phi0/Mp)
            dv = M**4*norm_deriv(phi0B/phi0,alpha,phi0/Mp)/phi0
            eps10 = norm_eps1(phi0B/phi0,alpha,phi0/Mp)

            mphiBoehm2 = 2 * M**4 / (phi0)**2 
            AB=mp.sqrt(80*alpha)*M**2/(phi0)
            lambdaB=Mp**3*mp.sqrt(alpha/5)*M**2/(phi0)**5

            return {'eps10':float(eps10), 'V0':float(V0), 'dv':float(dv)}
    
    with extradps(16):
        ret = aspic(lnRrad, ns, phi0, As)
    return ret

#print(tree(9e14, 2e-9, 0.9653, 0))
