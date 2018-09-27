"""
.. module:: functions
    :synopsis: some useful functions

"""

# define useful functions
import numpy as np
from scipy.special import erfc, spherical_jn, gamma
import mpmath
from mpmath import hyp2f1

def hyp2f1Reg(a,b,c,z):
    """
    return the regularized gauss hypergeometric function
    hypergeometric2F1Regularized in Mathematica
    simply equals hyp2f1(a,b,c,z)/Gamma(c)
    """
    return hyp2f1(a,b,c,z)/gamma(c)

def BesselJ(n, x):
    """
    return the spherical Bessel function of the first kind (J_n(x))
    """
    return sph_jn(n,x)[0][-1]

def xx(Mobs, M, sigma_lnM, BM=0):
    """
    return eq (16) of 1003.0841, the default systematic bias is set to 0.
    """
    return (np.log(Mobs)-BM-np.log(M))/np.sqrt(2)/sigma_lnM

def pp(Mobs, M, sigma_lnM):
    """
    return the probability of having an observed mass Mobs for which the true mass is M
    eqn 15 of 1003.0841
    """
    return np.exp(-xx(Mobs, M, sigma_lnM)**2.0)/np.sqrt(2.0*np.pi)/sigma_lnM

def selection_function(Mth, M, sigma_lnM):
    """
    return the selection function P_i(M)
    """
    return 0.5*erfc((np.log(Mth)-np.log(M))/(np.sqrt(2*sigma_lnM**2.0)))

def cubic_top_hat(L, kx, ky=0.0, kz=0.0):
    """return the Fourier transform of a cubic box top hat function
    """
    return np.power(L, 3.0)*np.sinc(kx*L/2.0)*np.sinc(ky*L/2.0)*np.sinc(kz*L/2.0)

def top_hat(k, R):
    """
    return the Fourier top hat
    """
    return (3.0*np.sin(k*R)-3.0*(k*R)*np.cos(k*R))/(k*R)**3.0

def dndlnM(M, rhom, sigma, dlnsigmainv_dlnM):
    """
    return the dndlnM value for the dark matter mass function fitting formula
    from Report of the DETF 2006
    """
    return -0.3*rhom*dlnsigmainv_dlnM*np.exp(-np.abs(np.log(1./sigma)+0.64)**3.82)/M

def Hermite3(x):
    """
    """
    return x**3.0-3*x

def atanxy(x, y, degrees=0):
    """ANGLE CCW FROM x-axis - from Dan Coe's fisher code"""
    theta = np.arctan(divsafe(y, x, inf=1e30, nan=0))
    theta = np.where(np.less(x, 0), theta + np.pi, theta)
    theta = np.where(np.logical_and(np.greater(x, 0), np.less(y, 0)), theta + 2*np.pi, theta)
    if degrees:
        theta = theta * 180. / np.pi
    return theta

def divsafe(a, b, inf=np.Inf, nan=np.NaN):
    """a / b with a / 0 = inf and 0 / 0 = nan"""
    a = np.array(a).astype(float)
    b = np.array(b).astype(float)
    asgn = np.greater_equal(a, 0) * 2 - 1.
    bsgn = np.greater_equal(b, 0) * 2 - 1.
    xsgn = asgn * bsgn
    sgn = np.where(b, xsgn, asgn)
    sgn = np.where(a, xsgn, bsgn)
    babs = np.clip(abs(b), 1e-200, 1e9999)
    bb = bsgn * babs
    #return where(b, a / bb, where(a, Inf, NaN))
    return np.where(b, a / bb, np.where(a, sgn*inf, nan))

# functions for returning fitting formulae that correct for baryonic effect, from 1402.4461
params=[[
         [-0.0881, 0.0881, -14.4100, 0.4280],
         [-0.0872, 0.0872, -13.6339, 0.3509]
         ],[
         [-0.0995, 0.0995, -14.7619, 0.3603],
         [-0.0993, 0.0995, -14.7619, 0.3603]
         ]]

def log10_fratio(M_DMonly, A, B, C, D):
    return A+B/(1.+np.exp(-(np.log10(M_DMonly)+C)/D))

def fratio(which, M_DMonly):
    """
    get the parameters using params[z][which]: currently which=0,1 are the only one supported referring to AGN8.5, M_200 and AGN8.0, M_200
    """
    prs=params[0][which]
    lvalue=log10_fratio(M_DMonly, prs[0], prs[1], prs[2], prs[3])
    return 10.0**lvalue

def tinker_f(sigma):
    """
    tinker f(sigma) Eqn (3), this is for Delta=200, z=0
    """
    A=0.186
    a=1.47
    b=2.57
    c=1.19
    return A*((sigma/b)**(-a)+1)*np.exp(-c/sigma**2.0)

def dlnsiginv(M):
    """
    expression for dlnsigma^-1/dlog10M using our previous simulation cosmology and fitting functions
    """
    return 6.51632+(269822.-14140.*np.log(M))/np.log(M)**3.0

def lnsigma(M):
    return 225.155+(58591.-6140.93*np.log(M))/(np.log(M)**2.0)-np.log(M**2.83)

def tinker_mf(M):
    rhomean=7.5e10; # in h^2M_sun/Mpc^3
    return tinker_f(np.exp(lnsigma(M)))*rhomean*dlnsiginv(M)/M

# end of baryonic effect in HMF from 1402.4461
def deg2rad(deg):
    return 0.01745329250*deg

def degsq2rad(degsq=41253, f=1.0):
    """
    convert degsq to radsq
    """
    return deg2rad(degsq)**2.0*f

def delta(i,j):
    if (i==j):
        return 1.
    else:
        return 0.
