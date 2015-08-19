#!/usr/bin/env python

import numpy as np
import inspect

RAD = np.pi/180.0
DEG = 180.0/np.pi
JD2000 = 2451545.00

def JCentury (JD):
    return (JD - JD2000)/36525.0
def T2JD(T):
    return 36525.0*T + JD2000

def eclipticObliquity (JD, debug=False):

    T = JCentury(JD)
    U = T / 100.0

    correction = - 4680.93 * U \
        -    1.55 * U**2\
        + 1999.25 * U**3\
        -   51.38 * U**4\
        -  249.67 * U**5\
        -   39.05 * U**6\
        +    7.12 * U**7\
        +   27.87 * U**8\
        +    5.79 * U**9\
        +    2.45 * U**10
    
    epsilon0 = 23.0 + 26.0/60.0 + (21.488 + correction)/3600.0
    if debug:
        print "%-24s %f" % (inspect.stack()[0][3], epsilon0)
    return epsilon0



def reduceAngle(angle):
    angle = angle - int(angle/360.0)*360
    if angle < 0.0:
        angle += 360.0
    return angle


# ------------------------------------------------------------------
# ------------------------------------------------------------------
def JD2epoch(JD):
    epoch = 2000 + (JD - JD2000)/365.25
    return epoch
def epoch2JD(epoch):
    jd = 365.25*(epoch - 2000.0) + JD2000
    return jd
def T2epoch(T):
    return JD2epoch(T2JD(T))

# --------------------------------------------------------------------
# --------------------------------------------------------------------
def sunGeoMeanLong (T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    L0 = reduceAngle(280.46646 + 36000.76983*T + 0.0003032*T**2.0)
    if debug:
        print "%-24s %f" % (inspect.stack()[0][3], L0)
    return L0


# --------------------------------------------------------------------
# --------------------------------------------------------------------
def sunMeanAnom (T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    M = reduceAngle(357.52911 + 35999.05029*T + 0.0001537*T**2.0)
    if debug:
        print "%-24s %f" % (inspect.stack()[0][3], M)
    return M



# --------------------------------------------------------------------
# --------------------------------------------------------------------
def sunEquationOfCenter (T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    M = sunMeanAnom (T)
    C = reduceAngle((1.914602 - 0.004817*T - 0.000014*T**2.0 )*np.sin(RAD*M)  \
                        + ( 0.019993 - 0.000101*T )*np.sin(RAD*2.0*M)  \
                        +  0.000289*np.sin(RAD*3.0*M))
    if debug:
        print "%-24s %f" % (inspect.stack()[0][3], C)
    return C




# --------------------------------------------------------------------
# --------------------------------------------------------------------
def sunTrueLongitude(T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    L0 = sunGeoMeanLong(T)
    C  = sunEquationOfCenter(T)
    ret = reduceAngle(L0 + C)
    if debug:
        print "%-24s %f" % (inspect.stack()[0][3], ret)
    
    return ret



# --------------------------------------------------------------------
# --------------------------------------------------------------------
def sunTrueLongJ2000(T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    epoch = T2epoch(T)
    trueLong = sunTrueLongitude(T)
    trueLongJ2000 = reduceAngle(trueLong - 0.01397*(epoch - 2000.0))
    if debug:
        print "%-24s %f" % (inspect.stack()[0][3], trueLongJ2000)
    return trueLongJ2000



# --------------------------------------------------------------------
# --------------------------------------------------------------------
def sunRAdec(T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    trueLongJ2000 = RAD * sunTrueLongJ2000(T)
    epsilon0      = RAD * eclipticObliquity(T)
    alpha = np.atan2(np.cos(epsilon0)*np.sin(trueLongJ2000), np.cos(trueLongJ2000))
    delta = np.asin( np.sin(epsilon0)*np.sin(trueLongJ2000) )
    alpha = reduceAngle(DEG*alpha)
    delta = DEG*delta
    if debug:
        print "%-24s %f %f" % (inspect.stack()[0][3], alpha, delta)
    return alpha, delta



# --------------------------------------------------------------------
# --------------------------------------------------------------------
def sunAppRAdec (T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    lamb = RAD * sunAppLongitude(T)
    Omega = RAD* (125.04 - 1934.136*T)
    epsilon0 = eclipticObliquity(T)
    epsilon = RAD * ( epsilon0 + 0.00256*np.cos(Omega) )
    alpha = np.atan2(np.cos(epsilon)*np.sin(lamb), np.cos(lamb))
    delta = np.asin( np.sin(epsilon)*np.sin(lamb) )
    alpha = reduceAngle(DEG*alpha)
    delta = DEG*delta
    if debug:
        print "%-24s %f %f" % (inspect.stack()[0][3], alpha, delta)
    return alpha, delta


def moonMeanLong(T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    Lp = reduceAngle(218.3164477 + 481267.88123421*T - 0.0015786*T*T + T*T*T/538841.0 - T**4/65194000.0)
    if debug:
        print "%-24s %f" % (inspect.stack()[0][3], Lp)
    return Lp

def moonMeanElong(T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    D  = reduceAngle(297.8501921 + 445267.1114034*T - 0.0018819*T*T + T*T*T/545868.0 - T**4/113065000.0)
    if debug:
        print "%-24s %f" % (inspect.stack()[0][3], D)
    return D

def moonMeanAnom(T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    Mp = reduceAngle(134.9633964 + 477198.8675055*T + 0.0087414*T*T + T*T*T/69699.0 + T**4/14712000.0)
    if debug:
        print "%-24s %f" % (inspect.stack()[0][3], Mp)
    return Mp

def moonArgOfLat(T, isJD=False, debug=False):
    if isJD:
        T = JCentury(T)
    F = reduceAngle(93.2720950 + 483202.0175233*T - 0.0036539*T*T - T*T*T/3526000.0 + T**4/863310000.0)
    if debug:
        print "%-24s %f" % (inspect.stack()[0][3], F)
    return F


def lunarPosition(T, isJD=False, debug=False):
    if isJD:
        T  = JCentury(T)
    Lp = RAD*moonMeanLong(T, debug=debug)
    D  = RAD*moonMeanElong(T, debug=debug)
    M  = RAD*sunMeanAnom(T, debug=debug)
    Mp = RAD*moonMeanAnom(T, debug=debug)
    F  = RAD*moonArgOfLat(T, debug=debug)
    A1 = RAD*reduceAngle((119.75 + 131.849*T))
    A2 = RAD*reduceAngle((53.09  + 479264.290*T))
    A3 = RAD*reduceAngle((313.45 + 481266.484*T))
    E  = RAD*reduceAngle(1.0 - 0.002516*T - 0.0000074*T*T)
    if debug:
        print "%-24s A1,A2,A3,E %f %f %f %f" % (inspect.stack()[0][3], DEG*A1, DEG*A2, DEG*A3, DEG*E)
    
    sigL = 3958.0*np.sin(A1) + 1962.0*np.sin(Lp - F) + 318.0*np.sin(A2)
    sigR = 0.0
    for arr in table47a():
        cD, cM, cMp, cF, sl, sr = arr
        if False: #cM in (1,-1):
            sl *= E
            sr *= E
        if False: #cM in (2, -2):
            sl *= E*E
            sr *= E*E
        arg = cD*D + cM*M + cMp*Mp + cF*F
        sigL += sl*np.sin(arg)
        sigR += sr*np.cos(arg)
    
    sigB = -2235.0*np.sin(Lp) \
        + 382.0*np.sin(A3) \
        + 175.0*np.sin(A1 - F) \
        + 175.0*np.sin(A1 + F) \
        + 127.0*np.sin(Lp - Mp) \
        - 115.0*np.sin(Lp + Mp)
    
    for arr in table47b():
        cD, cM, cMp, cF, sl = arr
        if False: #cM in (1,-1):
            print cM
            sl *= E
        if False: #cM in (2, -2):
            print cM
            sl *= E*E
        arg = cD*D + cM*M + cMp*Mp + cF*F
        sigB += sl*np.sin(arg)

    if debug:
        print "sL,B,R = ", -1127527, -3229126, -16590875
        print "%-24s sigL,B,R %f %f %f" % (inspect.stack()[0][3], sigL, sigB, sigR)
        
    lamb = reduceAngle(DEG*Lp + sigL/1000000.0)
    beta = sigB/1000000.0
    delta = 385000.56 + sigR/1000.0
    if debug:
        print "%-24s lamb,beta,delt %f %f %f" % (inspect.stack()[0][3], lamb, beta, delta)
    
    return lamb, beta, delta


def moonIllumFrac(T, isJD=False, debug=False):

    if isJD:
        T = JCentury(T)

    highPrecision = False
    
    if highPrecision:
        
        # high precision
        lamb, beta, earth_moon_dist = lunarPosition(T)
    
        alpha, delta     = sunAppRAdec(T)
        lamb0, _         = eq2ecl(RAD*alpha, RAD*delta)
        cos_psi = np.cos(beta)*np.cos(RAD*lamb - lamb0)
        psi = np.acos(cos_psi)
        Ro = 1.5e11
        tani = (Ro*np.sin(psi))/(earth_moon_dist - Ro*cos_psi)
        i = np.atan(tani)
    
    else:
        D  = RAD*moonMeanElong(T, debug=debug)
        M  = RAD*sunMeanAnom(T, debug=debug)
        Mp = RAD*moonMeanAnom(T, debug=debug)
        i  = reduceAngle(180.0
                         - DEG*D 
                         - 6.289*np.sin(Mp) 
                         + 2.100*np.sin(M) 
                         - 1.274*np.sin(2*D - Mp) 
                         - 0.658*np.sin(2*D) 
                         - 0.214*np.sin(2*Mp) 
                         - 0.110*np.sin(D))
    
    k = (1.0 + np.cos(RAD*i))/2.0
    if debug:
        print "%-24s %f %f" % (inspect.stack()[0][3], i, k)

    
    return k






# ------------------------------------------------------------------
# ------------------------------------------------------------------
def calendar2JD (Y, M, D, H=0, min=0, S=0):

    HpD = 24.0
    minpD = HpD*60.0
    SpD = minpD*60.0

    if ( M <= 2 ):
        Y -= 1
        M += 12
    
    A = int(Y/100)
    B =  2 - A + int(A/4) 

    (y,m,d) = (1582, 10, 4)
    if (Y<y or 
        (Y==y and M<m) or
        (Y==y and M==m and D<=4)):
        B = 0
    
    JD = int(365.25*(Y + 4716)) + int(30.6001*(M+1)) + D + B - 1524.5
    JD += H/HpD + min/minpD + S/SpD
    
    return JD






# ------------------------------------------------------------------
# ------------------------------------------------------------------
def JD2calendar (JD):

    JD += 0.5
    Z = int (JD)     # integer part
    F = JD - Z      # decimal part

    alpha = int( (Z - 1867216.25)/36524.25 )
    A = A if ( Z < 2299161 ) else Z + 1 + alpha - int(alpha/4)

    B = A + 1524
    C = int( (B - 122.1)/365.25 )
    D = int( 365.25*C )
    E = int( (B-D)/30.6001 )

    mday  = B - D - int(30.6001*E) + F
    mon   = E-1    if (E < 14) else (E-13)
    year  = C-4716 if (mon > 2) else (C-4715)

    hour = 24.0*F
    H = int(hour)
    min = (hour - H)*60.0
    Min = int(min)
    s = (min - Min)*60.0

    return (year, mon, mday, H, Min, s)



# ------------------------------------------------------------------
# ------------------------------------------------------------------
def calendar2epoch(Y, M, D, H=0, min=0, S=0):
    JD = calendar2JD(Y,M,D, H,min,S)
    epoch = 2000.0 + 100.0*JCentury(JD)
    return epoch



# ------------------------------------------------------------------
# ------------------------------------------------------------------
def epoch2calendar (epoch):
    jd = (epoch - 2000.0) * 365.25
    JD = JD2000 + jd
    (Y, M, D, H, min, S) = JD2calendar(JD)
    return (Y, M, D, H, min, S)



# ----------------------------------------------------------------
# ----------------------------------------------------------------
def yearDay(Y, M, D):
    is_leap = 1 if (not (Y % 4) and (Y % 400) ) else 0
    K = 1 if (is_leap) else 2
    yday = int(275.0*M/9.0) - K*int( (M+9.0)/12.0 ) + D - 30
    return yday




####################################################################
# 
#  12      12       12       12 
# 
###################################################################
def greenwichSidereal0hUT (JD):
    Y, M, D, H, m, S = JD2calendar(JD)
    JDmidnight = calendar2JD(Y, M, D, 0, 0, 0)
    T = JCentury(JDmidnight)
    theta0 = 100.46061837 + 36000.770053608*T + 0.0003879330*T**2 - T**3/38710000.0
    return reduceAngle(theta0)



def greenwichSidereal (JD):
    T = JCentury(JD)
    theta0 = 280.46061837 + 360.98564736629*(JD - JD2000) + 0.0003879330*T**2 - T**3/38710000.0
    return reduceAngle(theta0)























def ecl2eq(lamb, beta, JD):
    epsilon = RAD*eclipticObliquity(JD)
    numerator = np.sin(lamb)*np.cos(epsilon) - \
        np.tan(beta)*np.sin(epsilon)
    denominator = np.cos(lamb)

    alpha = np.atan2 (numerator, denominator)
    delta = np.asin( np.sin(beta)*np.cos(epsilon) +
                        np.cos(beta)*np.sin(epsilon)*np.sin(lamb) )

    return alpha, delta


def eq2ecl (alpha, delta, JD):
    epsilon = RAD*eclipticObliquity(JD)

    numerator = np.sin(alpha)*np.cos(epsilon) + np.tan(delta)*np.sin(epsilon)
    denominator = np.cos(alpha)

    lamb = np.atan2 (numerator, denominator)
    beta = np.asin( np.sin(delta)*np.cos(epsilon) -
                       np.cos(delta)*np.sin(epsilon)*sin(alpha) )

    return lamb, beta



#------------------------------------------------------------------
# function: mag2flux()
# Purpose: Get the flux from a star given its magnitude
# Req'd parameters: filter = photometric filter
#                   magnitude = self-expl.
#------------------------------------------------------------------
def mag2flux (filt, mag, exptime=1.0, radius=1.0/np.sqrt(np.pi)):

    area = np.pi*radius**2
    
    # get the flux
    # http://www.astro.utoronto.ca/~patton/astro/mags.html#flux
    #1 Jy = 10^-23 erg sec^-1 cm^-2 Hz^-1
    #1 Jy = 1.51e7 photons sec^-1 m^-2 (dlambda/lambda)^-1
    
    filter_specs = {
        'U' : [0.36, 0.15, 1810], 
        'B' : [0.44, 0.22, 4260], 
        'V' : [0.55, 0.16, 3640], 
        'R' : [0.64, 0.23, 3080], 
        'I' : [0.79, 0.19, 2550], 
        'J' : [1.26, 0.16, 1600], 
        'H' : [1.60, 0.23, 1080],
        'K' : [2.22, 0.23,  670],
        'g' : [0.52, 0.14, 3730],
        'r' : [0.67, 0.14, 4490],
        'i' : [0.79, 0.16, 4760],
        'z' : [0.91, 0.13, 4810]
        }

    if filt not in filter_specs:
        print "Warning Filter "+filt+" not in database."
        return 0.0
    
    # variable names  mag_flux_Jy  '_Jy'   --> mag_flux is *in* Janskys
    #                 photon_Flux_per_Jy   --> rate *per* Jansky
    lamb, dlambdaOverLambda, mag0_flux_Jy = filter_specs[filt]

    mag_flux_Jy = mag0_flux_Jy * 10**(-0.4*mag)

    photonFlux_per_Jy = 1.51e7 * dlambdaOverLambda

    mag_flux_phot = mag_flux_Jy * photonFlux_per_Jy * exptime * area

    return mag_flux_phot # photons per s per m^2 (if no exptime,area given)



def table47a():
    table = [
        # cD  cM  cMp  CF            sL             sR
        [0,   0,   1,   0,     6288774,     -20905355],
        [2,   0,  -1,   0,     1274027,      -3699111],
        [2,   0,   0,   0,      658314,      -2955968],
        [0,   0,   2,   0,      213618,       -569925],
        [0,   1,   0,   0,     -185116,         48888],
        [0,   0,   0,   2,     -114332,         -3149],
        [2,   0,  -2,   0,       58793,        246158],
        [2,  -1,  -1,   0,       57066,       -152138],
        [2,   0,   1,   0,       53322,       -170733],
        [2,  -1,   0,   0,       45758,       -204586],
        [0,   1,  -1,   0,      -40923,       -129620],
        [1,   0,   0,   0,      -34720,        108743],
        [0,   1,   1,   0,      -30383,        104755],
        [2,   0,   0,  -2,       15327,         10321],
        [0,   0,   1,   2,      -12528,             0],
        [0,   0,   1,  -2,       10980,         79661],
        [4,   0,  -1,   0,       10675,        -34782],
        [0,   0,   3,   0,       10034,        -23210],
        [4,   0,  -2,   0,       8548 ,        -21636],
        [2,   1,  -1,   0,      -7888 ,         24208],
        [2,   1,   0,   0,      -6766 ,         30824],
        [1,   0,  -1,   0,      -5163 ,         -8379],
        [1,   1,   0,   0,       4987 ,        -16675],
        [2,  -1,   1,   0,       4036 ,        -12831],
        [2,   0,   2,   0,       3994 ,        -10445],
        [4,   0,   0,   0,       3861 ,        -11650],
        [2,   0,  -3,   0,       3665 ,         14403],
        [0,   1,  -2,   0,      -2689 ,         -7003],
        [2,   0,  -1,   2,      -2602 ,             0],
        [2,  -1,  -2,   0,       2390 ,         10056],
        [1,   0,   1,   0,      -2348 ,          6322],
        [2,  -2,   0,   0,       2236 ,         -9884],
        
        [0,   1,   2,   0,      -2120 ,          5751],
        [0,   2,   0,   0,      -2069 ,             0],
        [2,  -2,  -1,   0,       2048 ,         -4950],
        [2,   0,   1,  -2,      -1773 ,          4130],
        [2,   0,   0,   2,      -1595 ,             0],
        [4,  -1,  -1,   0,       1215 ,         -3958],
        [0,   0,   2,   2,      -1110 ,             0],
        [3,   0,  -1,   0,       -892 ,          3258],
        [2,   1,   1,   0,       -810 ,          2616],
        [4,  -1,  -2,   0,        759 ,         -1897],
        [0,   2,  -1,   0,       -713 ,         -2117],
        [2,   2,  -1,   0,       -700 ,          2354],
        [2,   1,  -2,   0,        691 ,             0],
        [2,  -1,   0,  -2,        596 ,             0],
        [4,   0,   1,   0,        549 ,         -1423],
        [0,   0,   4,   0,        537 ,         -1117],
        [4,  -1,   0,   0,        520 ,         -1571],
        [1,   0,  -2,   0,       -487 ,         -1739],
        [2,   1,   0,  -2,       -399 ,             0],
        [0,   0,   2,  -2,       -381 ,         -4421],
        [1,   1,   1,   0,        351 ,             0],
        [3,   0,  -2,   0,       -340 ,             0],
        [4,   0,  -3,   0,        330 ,             0],
        [2,  -1,   2,   0,        327 ,             0],
        [0,   2,   1,   0,       -323 ,          1165],
        [1,   1,  -1,   0,        299 ,             0],
        [2,   0,   3,   0,        294 ,             0],
        [2,   0,  -1,  -2,          0 ,          8752],
        ]
    return table


def table47b():
    table = [
        [0,   0,  0,   1,    5128122],
        [0,   0,  1,   1,     280602],
        [0,   0,  1,  -1,     277693],
        [2,   0,  0,  -1,     173237],
        [2,   0, -1,   1,      55413],
        [2,   0, -1,  -1,      46271],
        [2,   0,  0,   1,      32573],
        [0,   0,  2,   1,      17198],
        [2,   0,  1,  -1,       9266],
        [0,   0,  2,  -1,       8822],
        [2,  -1,  0,  -1,       8216],
        [2,   0, -2,  -1,       4324],
        [2,   0,  1,   1,       4200],
        [2,   1,  0,  -1,      -3359],
        [2,  -1, -1,   1,       2463],
        [2,  -1,  0,   1,       2211],
        [2,  -1, -1,  -1,       2065],
        [0,   1, -1,  -1,      -1870],
        [4,   0, -1,  -1,       1828],
        [0,   1,  0,   1,      -1794],
        [0,   0,  0,   3,      -1749],
        [0,   1, -1,   1,      -1565],
        [1,   0,  0,   1,      -1491],
        [0,   1,  1,   1,      -1475],
        [0,   1,  1,  -1,      -1410],
        [0,   1,  0,  -1,      -1344],
        [1,   0,  0,  -1,      -1335],
        [0,   0,  3,   1,       1107],
        [4,   0,  0,  -1,       1021],
        [4,   0, -1,   1,        833],
        
        [0,   0,  1,  -3,        777],
        [4,   0, -2,   1,        671],
        [2,   0,  0,  -3,        607],
        [2,   0,  2,  -1,        596],
        [2,  -1,  1,  -1,        491],
        [1,   0, -2,   1,       -451],
        [0,   0,  3,  -1,        439],
        [2,   0,  2,   1,        422],
        [2,   0, -3,  -1,        421],
        [2,   1, -1,   1,       -366],
        [2,   1,  0,   1,       -351],
        [4,   0,  0,   1,        331],
        [2,  -1,  1,   1,        315],
        [2,  -2,  0,  -1,        302],
        [0,   0,  1,   3,       -283],
        [2,   1,  1,  -1,       -229],
        [1,   1,  0,  -1,        223],
        [1,   1,  0,   1,        223],
        [0,   1, -2,  -1,       -220],
        [2,   1, -1,  -1,       -220],
        [1,   0,  1,   1,       -185],
        [2,  -1, -2,  -1,        181],
        [0,   1,  2,   1,       -177],
        [4,   0, -2,  -1,        176],
        [4,  -1, -1,  -1,        166],
        [1,   0,  1,  -1,       -164],
        [4,   0,  1,  -1,        132],
        [1,   0, -1,  -1,       -119],
        [4,  -1,  0,  -1,        115],
        [2,  -2,  0,   1,        107],
        ]
    return table


if __name__ == '__main__':
    JD = 2448724.5
    lamb, beta, delt = lunarPosition(JD, debug=True)
    f = moonIllumFrac(JD, isJD=True, debug=True)

    cal = JD2calendar(JD)
    jd  = calendar2JD(*cal)
    print "Calendar: ", cal
    print "JD: ", jd

    epoch = JD2epoch(JD)
    jd    = epoch2JD(epoch)
    print "Epoch: ", epoch
    print "JD: ", jd
    
    print lamb, beta, delt
    print f
