from typing import NamedTuple

import numpy as np
import numpy.typing as npt

from leorbit.mathematics import Number
from leorbit.time import Timestamp, TimeInterval

### SGP4 ALGORITHM ###
### ############## ###

class OrbitalElementsComputeTuple(NamedTuple):
    n: int | float
    """Mean motion [rad/min]"""

    i: int | float
    """Inclination [rad]"""

    e: int | float
    """Eccentricity [1]"""

    argp: int | float
    """Argument of pericenter [rad]"""

    raan: int | float
    """Right ascension of ascending node [rad]"""

    M: int | float
    """Mean anomaly [rad]"""

    bstar: int | float
    """BSTAR drag term [1/earthRadii]"""

PI = np.float64(3.141592653589793)
TWOPI = np.float64(6.283185307179586)
HALFPI = np.float64(.5707963267948966)

XJ3 = -2.53881e-6
XKMPER = 6378.135 # Earth equatorial radius [km]
XMPER = XKMPER  * 1000
AE_MIN2M_S = XMPER / 60
GE = 398600.8 # Earth gravitational constant [km^3/s^2]
CK2 = 1.0826158e-3 / 2
CK4 = -3 * -1.65597e-6 / 8
S = 1 + 78 / XKMPER
QO = 1 + 120 / XKMPER
XKE = np.sqrt(3600 * GE / XKMPER**3)
QOMS2T = (QO - S)**4

class RawSGP4Output(NamedTuple):
    """Position and velocity vectors in GCRF frame, as returned by the SGP4 algorithm.
    Position coordinates are given in meters, velocity coordinates are given in meters per second."""

    x: npt.NDArray[np.float64] # [m]
    y: npt.NDArray[np.float64] # [m]
    z: npt.NDArray[np.float64] # [m]
    vx: npt.NDArray[np.float64] # [m/s]
    vy: npt.NDArray[np.float64] # [m/s]
    vz: npt.NDArray[np.float64] # [m/s]

def sgp4(
    elements_sat0: OrbitalElementsComputeTuple,
    tsince: npt.NDArray[np.float64] # [min]
) -> RawSGP4Output:
    """Python + Numpy implementation of SGP4 algorithm.
    
    Arguments:
    ----------
    - `sat0_n` Satellite's mean motion in given GP data. Units: `rad/min`
    - `sat0_i` Satellite's inclination in given GP data. Units: `rad`
    - `sat0_e` Satellite's eccentricity in given GP data. Units: `1`
    - `sat0_omega` Satellite's argument of periaster in given GP data. Units: `rad`
    - `sat0_node` Satellite's right ascension of ascending node in given GP data. Units: `rad`
    - `sat0_M` Satellite's mean anomaly in given GP data. Units `rad`
    - `bstar` Satellite's BSTAR coefficient in given GP data. Units `1/earthRadii`
    - `tsince` Time since GP data epoch. Units: `min`

    """
    sat0_n = np.float64(elements_sat0.n)
    sat0_i = np.float64(elements_sat0.i)
    sat0_e = np.float64(elements_sat0.e)
    sat0_omega = np.float64(elements_sat0.argp)
    sat0_node = np.float64(elements_sat0.raan)
    sat0_M = np.float64(elements_sat0.M)
    bstar = np.float64(elements_sat0.bstar)

    temp2 = XKE / sat0_n
    a1 = np.cbrt(temp2)**2
    cosio = np.cos(sat0_i)
    theta2 = cosio**2
    x3thm1 = 3*theta2 - 1
    eosq = sat0_e**2
    betao2 = 1 - eosq
    betao = np.sqrt(betao2)
    del1 = 1.5 * CK2 * x3thm1 / (a1**2 * betao * betao2)
    ao = a1 * ( 1 - del1*(1/3 + del1 * (1 + 134/81 * del1)))
    delo = 1.5 * CK2 * x3thm1 / (ao**2 * betao * betao2)
    xnodp = sat0_n / (1 + delo)
    aodp = ao / (1 - delo)

    isimp = 0
    if (aodp * (1 - sat0_e)) < (220.0 / XKMPER + 1):
        isimp = 1

    s4 = S
    qoms24 = QOMS2T
    perige = (aodp * (1 - sat0_e) - 1) * XKMPER
    if perige < 156.0:
        s4 = perige - 78.0
        if perige <= 98.0:
            s4 = 20.0
        qoms24 = ((120.0 - s4) / XKMPER)**4
        s4 = s4 / XKMPER + 1

    pinvsq = 1 / ( aodp**2 * betao2**2)
    tsi = 1 / (aodp - s4)
    eta = aodp * sat0_e * tsi
    etasq = eta**2
    eeta = sat0_e * eta
    psisq = np.fabs(1 - etasq)
    coef = qoms24 * tsi**4
    coef1 = coef / psisq**3.5
    c2 = coef1 * xnodp * (aodp * (1 + 1.5 * etasq + eeta * (4 + etasq)) + .75 * CK2 * tsi / psisq * x3thm1 * (8 + 3 * etasq * (8 + etasq)))
    c1 = bstar * c2
    sinio = np.sin(sat0_i)
    a3ovk2 = -XJ3 / CK2
    c3 = coef * tsi * a3ovk2 * xnodp * sinio / sat0_e
    x1mth2 = 1 - theta2
    c4 = 2 * xnodp * coef1 * aodp * betao2 * (eta * (2 + .5 * etasq) + sat0_e * (.5 + 2 * etasq) - 2 * CK2 * tsi / (aodp * psisq) * ( -3 * x3thm1 * ( 1 - 2 * eeta + etasq * (1.5 - .5 * eeta)) + .75 * x1mth2 * (2 * etasq - eeta * (1 + etasq)) * np.cos(2 * sat0_omega)))
    c5 = 2 * coef1 * aodp * betao2 * (1 + 2.75 * (etasq + eeta) + eeta * etasq)
    theta4 = theta2**2
    temp1 = 3 * CK2 * pinvsq * xnodp
    temp2 = temp1 * CK2 * pinvsq
    temp3 = 1.25 * CK4 * pinvsq * pinvsq * xnodp
    xmdot = xnodp + .5 * temp1 * betao * x3thm1 + .0625 * temp2 * betao * (13 - 78 * theta2 + 137 * theta4)
    x1m5th = 1.0 - 5.0 * theta2
    omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4)
    xhdot1 = -temp1 * cosio
    xnodot = xhdot1 + (.5 * temp2 * (4 - 19 * theta2) + 2 * temp3 * (3 - 7 * theta2)) * cosio
    omgcof = bstar * c3 * np.cos(sat0_omega)
    xmcof = -2/3 * coef * bstar / eeta
    xnodcf = 3.5 * betao2 * xhdot1 * c1
    t2cof = 1.5 * c1
    xlcof = .125 * a3ovk2 * sinio * (3 + 5 * cosio) / (1 + cosio)
    aycof = .25 * a3ovk2 * sinio
    delmo = (1 + eta * np.cos(sat0_M))**3
    sinmo = np.sin(sat0_M)
    x7thm1 = 7 * theta2 - 1

    c1sq =  0
    d2 = 0
    temp =  0
    d3 = 0
    d4 = 0
    t3cof = 0
    t4cof = 0
    t5cof = 0
    if isimp == 0:
        c1sq = c1**2
        d2 = 4 * aodp * tsi * c1sq
        temp = d2 * tsi * c1 / 3
        d3 = (17 * aodp + s4) * temp
        d4 = .5 * temp * aodp * tsi * (221 * aodp + 31 * s4) * c1
        t3cof = d2 + 2 * c1sq
        t4cof = .25 * (3 * d3 + c1 * (12 * d2 + 10 * c1sq))
        t5cof = .2 * (3 * d4 + 12 * c1 * d3 + 6 * d2 * d2 + 15 * c1sq * (2 * d2 + c1sq))

    xmdf = sat0_M + xmdot * tsince
    omgadf = sat0_omega + omgdot * tsince
    xnoddf = sat0_node + xnodot * tsince
    omega = omgadf
    xmp = xmdf
    tsq = tsince**2
    xnode = xnoddf + xnodcf * tsq
    tempa = 1 - c1 * tsince
    tempe = bstar * c4 * tsince
    templ = t2cof * tsq

    delomg = 0
    delm = 0
    tcube = 0
    tfour = 0
    if isimp == 0:
        delomg = omgcof * tsince
        delm = xmcof*(((1 + eta * np.cos(xmdf))**3) - delmo)
        temp = delomg + delm
        xmp = xmdf + temp
        omega = omgadf - temp
        tcube = tsq * tsince
        tfour = tsince * tcube
        tempa = tempa - d2 * tsq - d3 * tcube - d4 * tfour
        tempe = tempe + bstar * c5 * (np.sin(xmp) - sinmo)
        templ = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof)

    a = aodp * tempa**2
    e = sat0_e - tempe
    xl = xmp + omega + xnode + xnodp * templ
    beta = np.sqrt(1 - e**2)
    xn = XKE / a**1.5

    axn = e * np.cos(omega)
    temp = 1 / (a * beta**2)
    xll = temp * xlcof * axn
    aynl = temp * aycof
    xlt = xl + xll
    ayn = e * np.sin(omega) + aynl

    capu = (xlt - xnode) % TWOPI
    temp2 = capu

    sinepw = 0
    cosepw = 0
    temp4 = 0
    temp5 = 0
    temp6 = 0
    epw = 0
    temp7 = 0
    for _ in range(10):
        sinepw = np.sin(temp2)
        cosepw = np.cos(temp2)
        temp3 = axn * sinepw
        temp4 = ayn * cosepw
        temp5 = axn * cosepw
        temp6 = ayn * sinepw
        epw = (capu - temp4 + temp3 - temp2) / (1 - temp5 - temp6) + temp2
        temp7 = temp2
        temp2 = epw
        cond_break = np.fabs(epw - temp7) <= 1e-6
        if np.all(cond_break):
            break

    ecose = temp5 + temp6
    esine = temp3 - temp4
    elsq = axn**2 + ayn**2
    temp = 1 - elsq
    pl = a * temp
    r = a * (1 - ecose)
    temp1 = 1 / r
    rdot = XKE * np.sqrt(a) * esine * temp1
    rfdot = XKE * np.sqrt(pl) * temp1
    temp2 = a * temp1
    betal = np.sqrt(temp)
    temp3 = 1 / (1 + betal)
    cosu = temp2 * (cosepw - axn + ayn * esine * temp3)
    sinu = temp2 * (sinepw - ayn - axn * esine * temp3)
    u = np.atan2(sinu, cosu)
    sin2u = 2 * sinu * cosu
    cos2u = 2 * cosu**2 - 1
    temp = 1 / pl
    temp1 = CK2 * temp
    temp2 = temp1 * temp

    rk = r * (1 - 1.5 * temp2 * betal * x3thm1) + .5 * temp1 * x1mth2 * cos2u # [radiiEarth]
    uk = u - .25 * temp2 * x7thm1 * sin2u
    xnodek = xnode + 1.5 * temp2 * cosio * sin2u
    xinck = sat0_i + 1.5 * temp2 * cosio * sinio * cos2u
    rdotk = rdot - xn * temp1 * x1mth2 * sin2u # [radiiEarth/min]
    rfdotk = rfdot + xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1) # [radiiEarth/min]
    
    c_raan = np.cos(xnodek)
    s_raan = np.sin(xnodek)
    c_i = np.cos(xinck)
    s_i = np.sin(xinck)
    c_theta = np.cos(uk)
    s_theta = np.sin(uk)

    u_x = c_raan * c_theta - s_raan * c_i * s_theta
    u_y = s_raan * c_theta + c_raan * c_i * s_theta
    u_z = s_i * s_theta

    v_x = - c_raan * s_theta - s_raan * c_i * c_theta
    v_y = - s_raan * s_theta + c_raan * c_i * c_theta
    v_z = s_i * c_theta

    gcrf_x = (rk * u_x) * XMPER
    gcrf_y = (rk * u_y) * XMPER
    gcrf_z = (rk * u_z) * XMPER

    gcrf_vel_x = (rdotk * u_x + rfdotk * v_x) * AE_MIN2M_S
    gcrf_vel_y = (rdotk * u_y + rfdotk * v_y) * AE_MIN2M_S
    gcrf_vel_z = (rdotk * u_z + rfdotk * v_z) * AE_MIN2M_S

    return RawSGP4Output(gcrf_x, gcrf_y, gcrf_z, gcrf_vel_x, gcrf_vel_y, gcrf_vel_z)