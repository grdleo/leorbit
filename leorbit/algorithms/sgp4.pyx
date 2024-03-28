from math import sin, cos, atan, sqrt, fabs, cbrt, atan2

import cython

cdef double PI = 3.141592653589793
cdef double TWOPI = 6.283185307179586
cdef double HALFPI = 1.5707963267948966

cdef double XJ3 = -2.53881e-6
cdef double XKMPER = 6378.135 # Earth equatorial radius [km]
cdef double XMPER = XKMPER  * 1000
cdef double AE_MIN2M_S = XMPER / 60
cdef double GE = 398600.8 # Earth gravitational constant [km^3/s^2]
cdef double CK2 = 1.0826158e-3 / 2
cdef double CK4 = -3 * -1.65597e-6 / 8
cdef double S = 1 + 78 / XKMPER
cdef double QO = 1 + 120 / XKMPER
cdef double XKE = sqrt(3600 * GE / XKMPER**3)
cdef double QOMS2T = (QO - S)**4

cpdef tuple algo(
    double sat0_n, # [rad/min]
    double sat0_i, # [rad]
    double sat0_e, # [1]
    double sat0_omega, # [rad]
    double sat0_node, # [rad]
    double sat0_M, # [rad]
    double bstar, # [1/earthRadii]
    double tsince # [min]
):
    
    cdef double temp2 = XKE / sat0_n
    cdef double a1 = cbrt(temp2)**2
    cdef double cosio = cos(sat0_i)
    cdef double theta2 = cosio**2
    cdef double x3thm1 = 3*theta2 - 1
    cdef double eosq = sat0_e**2
    cdef double betao2 = 1 - eosq
    cdef double betao = sqrt(betao2)
    cdef double del1 = 1.5 * CK2 * x3thm1 / (a1**2 * betao * betao2)
    cdef double ao = a1 * ( 1 - del1*(1/3 + del1 * (1 + 134/81 * del1)))
    cdef double delo = 1.5 * CK2 * x3thm1 / (ao**2 * betao * betao2)
    cdef double xnodp = sat0_n / (1 + delo)
    cdef double aodp = ao / (1 - delo)

    # Initialization
    # For perigee less than 220 kilometers, the isimp flag is set and
    # the equations are truncated to linear variation in sqrt a and
    # quadratic variation in mean anomaly.  Also, the c3 term, the
    # delta omega term, and the delta m term are dropped.

    cdef int isimp = 0
    if (aodp * (1 - sat0_e)) < (220.0 / XKMPER + 1):
        isimp = 1

    # For perigee below 156 km, the values of s and QOMS2T are altered.
    cdef double s4 = S
    cdef double qoms24 = QOMS2T
    cdef double perige = (aodp * (1 - sat0_e) - 1) * XKMPER
    if perige < 156.0:
        s4 = perige - 78.0
        if perige <= 98.0:
            s4 = 20.0
        qoms24 = ((120.0 - s4) / XKMPER)**4
        s4 = s4 / XKMPER + 1

    cdef double pinvsq = 1 / ( aodp**2 * betao2**2)
    cdef double tsi = 1 / (aodp - s4)
    cdef double eta = aodp * sat0_e * tsi
    cdef double etasq = eta**2
    cdef double eeta = sat0_e * eta
    cdef double psisq = fabs(1 - etasq)
    cdef double coef = qoms24 * tsi**4
    cdef double coef1 = coef / psisq**3.5
    cdef double c2 = coef1 * xnodp * (aodp * (1 + 1.5 * etasq + eeta * (4 + etasq)) + .75 * CK2 * tsi / psisq * x3thm1 * (8 + 3 * etasq * (8 + etasq)))
    cdef double c1 = bstar * c2
    cdef double sinio = sin(sat0_i)
    cdef double a3ovk2 = -XJ3 / CK2
    cdef double c3 = coef * tsi * a3ovk2 * xnodp * sinio / sat0_e
    cdef double x1mth2 = 1 - theta2
    cdef double c4 = 2 * xnodp * coef1 * aodp * betao2 * (eta * (2 + .5 * etasq) + sat0_e * (.5 + 2 * etasq) - 2 * CK2 * tsi / (aodp * psisq) * ( -3 * x3thm1 * ( 1 - 2 * eeta + etasq * (1.5 - .5 * eeta)) + .75 * x1mth2 * (2 * etasq - eeta * (1 + etasq)) * cos(2 * sat0_omega)))
    cdef double c5 = 2 * coef1 * aodp * betao2 * (1 + 2.75 * (etasq + eeta) + eeta * etasq)
    cdef double theta4 = theta2**2
    cdef double temp1 = 3 * CK2 * pinvsq * xnodp
    temp2 = temp1 * CK2 * pinvsq
    cdef double temp3 = 1.25 * CK4 * pinvsq * pinvsq * xnodp
    cdef double xmdot = xnodp + .5 * temp1 * betao * x3thm1 + .0625 * temp2 * betao * (13 - 78 * theta2 + 137 * theta4)
    cdef double x1m5th = 1.0 - 5.0 * theta2
    cdef double omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4)
    cdef double xhdot1 = -temp1 * cosio
    cdef double xnodot = xhdot1 + (.5 * temp2 * (4 - 19 * theta2) + 2 * temp3 * (3 - 7 * theta2)) * cosio
    cdef double omgcof = bstar * c3 * cos(sat0_omega)
    cdef double xmcof = -2/3 * coef * bstar / eeta
    cdef double xnodcf = 3.5 * betao2 * xhdot1 * c1
    cdef double t2cof = 1.5 * c1
    cdef double xlcof = .125 * a3ovk2 * sinio * (3 + 5 * cosio) / (1 + cosio)
    cdef double aycof = .25 * a3ovk2 * sinio
    cdef double delmo = (1 + eta * cos(sat0_M))**3
    cdef double sinmo = sin(sat0_M)
    cdef double x7thm1 = 7 * theta2 - 1

    cdef double c1sq =  0
    cdef double d2 = 0
    cdef double temp =  0
    cdef double d3 = 0
    cdef double d4 = 0
    cdef double t3cof = 0
    cdef double t4cof = 0
    cdef double t5cof = 0
    if isimp == 0:
        c1sq = c1**2
        d2 = 4 * aodp * tsi * c1sq
        temp = d2 * tsi * c1 / 3
        d3 = (17 * aodp + s4) * temp
        d4 = .5 * temp * aodp * tsi * (221 * aodp + 31 * s4) * c1
        t3cof = d2 + 2 * c1sq
        t4cof = .25 * (3 * d3 + c1 * (12 * d2 + 10 * c1sq))
        t5cof = .2 * (3 * d4 + 12 * c1 * d3 + 6 * d2 * d2 + 15 * c1sq * (2 * d2 + c1sq))

    cdef double xmdf = sat0_M + xmdot * tsince
    cdef double omgadf = sat0_omega + omgdot * tsince
    cdef double xnoddf = sat0_node + xnodot * tsince
    cdef double omega = omgadf
    cdef double xmp = xmdf
    cdef double tsq = tsince**2
    cdef double xnode = xnoddf + xnodcf * tsq
    cdef double tempa = 1 - c1 * tsince
    cdef double tempe = bstar * c4 * tsince
    cdef double templ = t2cof * tsq

    cdef double delomg = 0
    cdef double delm = 0
    cdef double tcube = 0
    cdef double tfour = 0
    if isimp == 0:
        delomg = omgcof * tsince
        delm = xmcof*(((1 + eta * cos(xmdf))**3) - delmo)
        temp = delomg + delm
        xmp = xmdf + temp
        omega = omgadf - temp
        tcube = tsq * tsince
        tfour = tsince * tcube
        tempa = tempa - d2 * tsq - d3 * tcube - d4 * tfour
        tempe = tempe + bstar * c5 * (sin(xmp) - sinmo)
        templ = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof)

    cdef double a = aodp * tempa**2
    cdef double e = sat0_e - tempe
    cdef double xl = xmp + omega + xnode + xnodp * templ
    cdef double beta = sqrt(1 - e**2)
    cdef double xn = XKE / a**1.5

    cdef double axn = e * cos(omega)
    temp = 1 / (a * beta**2)
    cdef double xll = temp * xlcof * axn
    cdef double aynl = temp * aycof
    cdef double xlt = xl + xll
    cdef double ayn = e * sin(omega) + aynl

    cdef double capu = (xlt - xnode) % TWOPI
    temp2 = capu

    cdef double sinepw = 0
    cdef double cosepw = 0
    cdef double temp4 = 0
    cdef double temp5 = 0
    cdef double temp6 = 0
    cdef double epw = 0
    cdef double temp7 = 0
    for i in range(10):
        sinepw = sin(temp2)
        cosepw = cos(temp2)
        temp3 = axn * sinepw
        temp4 = ayn * cosepw
        temp5 = axn * cosepw
        temp6 = ayn * sinepw
        epw = (capu - temp4 + temp3 - temp2) / (1 - temp5 - temp6) + temp2
        temp7 = temp2
        temp2 = epw
        if fabs(epw - temp7) <= 1e-6:
            break

    cdef double ecose = temp5 + temp6
    cdef double esine = temp3 - temp4
    cdef double elsq = axn**2 + ayn**2
    temp = 1 - elsq
    cdef double pl = a * temp
    cdef double r = a * (1 - ecose)
    temp1 = 1 / r
    cdef double rdot = XKE * sqrt(a) * esine * temp1
    cdef double rfdot = XKE * sqrt(pl) * temp1
    temp2 = a * temp1
    cdef double betal = sqrt(temp)
    temp3 = 1 / (1 + betal)
    cdef double cosu = temp2 * (cosepw - axn + ayn * esine * temp3)
    cdef double sinu = temp2 * (sinepw - ayn - axn * esine * temp3)
    cdef double u = atan2(sinu, cosu)
    cdef double sin2u = 2 * sinu * cosu
    cdef double cos2u = 2 * cosu**2 - 1
    temp = 1 / pl
    temp1 = CK2 * temp
    temp2 = temp1 * temp

    cdef double rk = r * (1 - 1.5 * temp2 * betal * x3thm1) + .5 * temp1 * x1mth2 * cos2u # [radiiEarth]
    cdef double uk = u - .25 * temp2 * x7thm1 * sin2u
    cdef double xnodek = xnode + 1.5 * temp2 * cosio * sin2u
    cdef double xinck = sat0_i + 1.5 * temp2 * cosio * sinio * cos2u
    cdef double rdotk = rdot - xn * temp1 * x1mth2 * sin2u # [radiiEarth/min]
    cdef double rfdotk = rfdot + xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1) # [radiiEarth/min]
    
    cdef double c_raan = cos(xnodek)
    cdef double s_raan = sin(xnodek)
    cdef double c_i = cos(xinck)
    cdef double s_i = sin(xinck)
    cdef double c_theta = cos(uk)
    cdef double s_theta = sin(uk)

    cdef double u_x = c_raan * c_theta - s_raan * c_i * s_theta,
    cdef double u_y = s_raan * c_theta + c_raan * c_i * s_theta,
    cdef double u_z = s_i * s_theta

    cdef double v_x = - c_raan * s_theta - s_raan * c_i * c_theta,
    cdef double v_y = - s_raan * s_theta + c_raan * c_i * c_theta,
    cdef double v_z = s_i * c_theta

    cdef double gcrf_x = (rk * u_x) * XMPER
    cdef double gcrf_y = (rk * u_y) * XMPER
    cdef double gcrf_z = (rk * u_z) * XMPER

    cdef double gcrf_vel_x = (rdotk * u_x + rfdotk * v_x) * AE_MIN2M_S
    cdef double gcrf_vel_y = (rdotk * u_y + rfdotk * v_y) * AE_MIN2M_S
    cdef double gcrf_vel_z = (rdotk * u_z + rfdotk * v_z) * AE_MIN2M_S

    return gcrf_x, gcrf_y, gcrf_z, gcrf_vel_x, gcrf_vel_y, gcrf_vel_z