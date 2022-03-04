import numpy as np
import math
import scipy.optimize as opt


class EOSparams:
    """
    h                       enthalpy of formation
    s                       entropy
    v0                      volume at 1bar, 298K
    n                       atoms per formula unit
    a0                      coefficient of thermal expansion
    K                       bulk modulus at 1bar, 298K
    dKdP                    first derivative of K
    dKdP2                   second derivative of K
    cp_a, cp_b, cp_c, cp_d  coefficients of heat capacity polynomial: a + bT + cT**-2 + dT**-(1/2)
    smax, vmax, Tc0         Landau theory parameters for phase transitions

    Most data from Holland & Powell (2011), some updated from somewhere else by Olivier
    """

    fayalite = {
        "h": -1477.510,
        "s": 151.0,
        "v0": 4.631,
        "n": 7,
        "a0": 2.82e-5,
        "K": 1256,
        "dKdP": 4.68,
        "dKdP2": -3.7e-3,
        "cp_a": 2.011e-1,
        "cp_b": 1.733e-5,
        "cp_c": -1960.6,
        "cp_d": -9.009e-1,
    }

    ringwoodite = {
        "h": -1477.510,
        "s": 140.0,
        "v0": 4.203,
        "n": 7,
        "a0": 2.22e-5,
        "K": 1977,
        "dKdP": 4.92,
        "dKdP2": -2.5e-3,
        "cp_a": 1.668e-1,
        "cp_b": 4.2610e-5,
        "cp_c": -1705.4,
        "cp_d": -5.414e-1,
    }

    quartz = {
        "h": -910.710,
        "s": 41.43,
        "v0": 2.269,
        "n": 3,
        "a0": 0,
        "K": 730,
        "dKdP": 6,
        "dKdP2": -8.2e-3,
        "smax": 4.95 / 1e3,
        "vmax": 1.188e-1,
        "Tc0": 847,
        "cp_a": 9.29e-2,
        "cp_b": -6.42e-7,
        "cp_c": -714.9,
        "cp_d": -0.7161,
    }

    coesite = {
        "h": -906.990,
        "s": 39.60,
        "v0": 2.064,
        "n": 3,
        "a0": 1.23e-5,
        "K": 979,
        "dKdP": 4.19,
        "dKdP2": -4.3e-3,
        "cp_a": 1.078e-1,
        "cp_b": -3.279e-6,
        "cp_c": -190.3,
        "cp_d": -1.0416,
    }

    stishovite = {
        "h": -876.720,
        "s": 24.0,
        "v0": 1.401,
        "n": 3,
        "a0": 1.58e-5,
        "K": 3090,
        "dKdP": 4.6,
        "dKdP2": -1.50e-3,
        "cp_a": 6.81e-2,
        "cp_b": 6.010e-6,
        "cp_c": -1978.2,
        "cp_d": -8.21e-2,
    }

    magnetite = {
        "h": -1114.510,
        "s": 146.9,
        "v0": 4.452,
        "n": 7,
        "a0": 3.71e-5,
        "K": 1857,
        "dKdP": 4.05,
        "dKdP2": -2.2e-3,
        "smax": 35.0,
        "vmax": 0.0,
        "Tc0": 848,
    }


def teosth(phase, pkbar, t, tref=298.15, **kwargs):
    """Tait equation of state"""

    params = ["s", "v0", "n", "a0", "K", "dKdP", "dKdP2"]
    s, v0, n, a0, K, dKdP, dKdP2 = [getattr(EOSparams, phase)[i] for i in params]

    theta = 10636 / (s / n + 6.44)
    u0 = theta / tref
    ksi0 = math.pow(u0, 2) * math.exp(u0) / math.pow((math.exp(u0) - 1), 2)
    a = (1 + dKdP) / (1 + dKdP + K * dKdP2)
    b = dKdP / K - dKdP2 / (1 + dKdP)
    c = (1 + dKdP + K * dKdP2) / (math.pow(dKdP, 2) + dKdP - K * dKdP2)
    u = theta / t
    pth = a0 * K * theta / ksi0 * (1 / (math.exp(u) - 1) - 1 / (math.exp(u0) - 1))
    intVdP = (
        pkbar
        * v0
        * (
            1
            - a
            + a
            * (
                np.sign((1 - b * pth)) * math.pow(abs(1 - b * pth), (1 - c))
                - np.sign((1 + b * (pkbar - pth)))
                * math.pow(abs((1 + b * (pkbar - pth))), (1 - c))
            )
            / (b * (c - 1) * pkbar)
        )
    )

    return [pth, a, b, c, intVdP]


def teosthV(v, pth, a, b, c, v0, pkbar, **kwargs):
    """Tait equation of state"""

    vteosth = v - v0 * (1 - a * (1 - (1 + b * (pkbar - pth)) ** -c))

    return vteosth


def landauhppressure(phase, pkbar, t):

    smax, vmax, tc0 = [getattr(EOSparams, phase)[i] for i in ["smax", "vmax", "Tc0"]]

    q02 = np.sqrt(1 - 298.15 / tc0)
    tc = tc0 + pkbar * vmax / smax
    if (tc - t) > 0:
        q2 = np.sqrt((tc - t) / tc0)
    else:
        q2 = 0
    gdistot = (
        smax
        * (tc0 * (q02 - 1 / 3 * q02 ** 3. + 1 / 3 * q2 ** 3.) - tc * q2 - t * (q02 - q2))
        + pkbar * vmax * q02
    )
    tc = tc0 + pkbar * 0 / smax
    if (tc - t) > 0:
        q2 = np.sqrt((tc - t) / tc0)
    else:
        q2 = 0
    gdisprind = smax * (
        tc0 * (q02 - 1 / 3 * q02 ** 3. + 1 / 3 * q2 ** 3.) - tc * q2 - t * (q02 - q2)
    )
    gdisprdep = gdistot - gdisprind

    return gdisprdep


def CpHoll(phase, t, tref=298.15):
    """Heat capacity from Holland"""

    a, b, c, d = [
        getattr(EOSparams, phase)[i] for i in ["cp_a", "cp_b", "cp_c", "cp_d"]
    ]

    CpHolland = (a * t + 0.5 * b * (t ** 2.) - c * (t ** -1.) + 2 * d * t ** 0.5) - (
        a * tref + 0.5 * b * (tref ** 2.) - c * (tref ** -1.) + 2 * d * tref ** 0.5
    )

    return CpHolland


def CpTHoll(phase, t, tref=298.15):
    """Heat capacity from Holland"""

    a, b, c, d = [
        getattr(EOSparams, phase)[i] for i in ["cp_a", "cp_b", "cp_c", "cp_d"]
    ]

    CpTHolland = (a * np.log(t) + b * t - 0.5 * c * t ** (-2.) - 2 * d * t ** (-0.5)) - (
        a * np.log(tref) + b * tref - 0.5 * c * tref ** (-2.) - 2 * d * tref ** (-0.5)
    )

    return CpTHolland


def landauhp1bar(phase, pkbar, t):

    smax, vmax, tc0 = [getattr(EOSparams, phase)[i] for i in ["smax", "vmax", "Tc0"]]

    q02 = np.sqrt(1 - 298.15 / tc0)
    tc = tc0 + pkbar * vmax / smax
    if (tc - t) > 0:
        q2 = np.sqrt((tc - t) / tc0)
    else:
        q2 = 0
    gdistot = (
        smax
        * (tc0 * (q02 - 1 / 3 * q02 ** 3. + 1 / 3 * q2 ** 3.) - tc * q2 - t * (q02 - q2))
        + pkbar * vmax * q02
    )

    return gdistot
