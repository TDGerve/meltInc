from . import eos
from .. import chemcalc as cc
import scipy.optimize as opt
import warnings as w
import numpy as np
import pandas as pd
from scipy.constants import R


def QFM_pressure(T_K, Pbar):

    pgpa = Pbar / 1e4

    # Pressure of transition SiO2 polymorphs
    pqzcoe = lambda pkbar: phaseTransition(
        pkbar, t=T_K, phase_1="quartz", phase_2="coesite"
    )
    pqzcoe = opt.fsolve(pqzcoe, 8)

    pcoestish = lambda pkbar: phaseTransition(
        pkbar, t=T_K, phase_1="coesite", phase_2="stishovite"
    )
    pcoestish = opt.fsolve(pcoestish, 8)

    # Pressure of transition Fe2SiO4 polymorphs
    pfaringHol_P = lambda pkbar: phaseTransition(
        pkbar, t=T_K, phase_1="fayalite", phase_2="ringwoodite"
    )
    pfaringHol_P = opt.fsolve(pfaringHol_P, 8)

    # VdP SiO2 polymorphs
    *_, vdpqz = eos.teosth(phase="quartz", pkbar=min(pgpa * 10, pqzcoe), t=T_K)
    qlandau = eos.landauhppressure(phase="quartz", pkbar=min(pgpa * 10, pqzcoe), t=T_K)
    dvdpSiO2Hol_P = vdpqz + qlandau
    if pgpa > 0.1 * pqzcoe:
        vdcoe = (
            eos.teosth(phase="coesite", pkbar=min(pgpa * 10, pcoestish), t=T_K)[4]
            - eos.teosth(phase="coesite", pkbar=pqzcoe, t=T_K)[4]
        )
        dvdpSiO2Hol_P = dvdpSiO2Hol_P + vdcoe
        if pgpa > 0.1 * pcoestish:
            vdpstish = (
                eos.teosth(phase="stishovite", pkbar=pgpa * 10, t=T_K)[4]
                - eos.teosth(phase="stishovite", pkbar=pcoestish, t=T_K)[4]
            )
            dvdpSiO2Hol_P = dvdpSiO2Hol_P + vdpstish

    # VdP Fe2SiO4 polymorphs
    *_, vdpfa = eos.teosth(phase="fayalite", pkbar=min(pgpa * 10, pfaringHol_P), t=T_K)
    dvdpFe2SiO4Hol_P = vdpfa
    if pgpa > 0.1 * pfaringHol_P:
        vdpring = (
            eos.teosth(phase="ringwoodite", pkbar=pgpa * 10, t=T_K)[4]
            - eos.teosth(phase="ringwoodite", pkbar=pfaringHol_P, t=T_K)[4]
        )
        dvdpFe2SiO4Hol_P = dvdpFe2SiO4Hol_P + vdpring

    # VdP magnetite
    vdpmt = eos.teosth(phase="magnetite", pkbar=pgpa * 10, t=T_K)[4]

    return 3e3 * dvdpSiO2Hol_P + 2e3 * vdpmt - 3e3 * dvdpFe2SiO4Hol_P


def phaseTransition(pkbar, t, phase_1, phase_2):

    results = []

    for phase in [phase_1, phase_2]:

        h = getattr(eos.EOSparams, phase)["h"]
        s = getattr(eos.EOSparams, phase)["s"] / 1e3

        Gibbs = (
            h + eos.CpHoll(phase=phase, t=t) - t * (s + eos.CpTHoll(phase=phase, t=t))
        )
        *_, vdp = eos.teosth(phase=phase, pkbar=pkbar, t=t)

        Gibbs = Gibbs + vdp

        if phase in ["quartz", "magnetite"]:

            Gibbs = Gibbs + eos.landauhp1bar(phase=phase, pkbar=pkbar, t=t)

        results.append(Gibbs)

    return results[0] - results[1]


def QFM_1bar(T_K):
    """
    calculate chemical potential of oxygen (fO2) at QFM a 1 bar. Equation from O'Neill 1987

    Parameters
    ----------
    T_K     list-like, float
        Temperature in Kelvin

    Returns
    -------
    Chemical potential (or Gibbs free energy per mole)
    """

    # if T_K is an integer
    try:
        if T_K < 900:
            w.warn("Temperature below 900K")
        if T_K > 1420:
            w.warn("Temperature above 1420K")
    except:
        pass
    # if T_K is list-like
    try:
        if (np.array(T_K) < 900).any():
            w.warn("Temperatures below 900K present")
        if (np.array(T_K) > 1420).any():
            w.warn("Temperatures above 1420K present")
    except:
        pass

    return -587474 + 1584.427 * T_K - 203.3164 * T_K * np.log(T_K) + 0.092710 * T_K ** 2


def fO2QFM(logshift, T_K, Pbar):
    """
    calculate oxygen fugacity at QFM offset by arbitraty log units.

    Parameters
    ----------
    log_units   int
        Log units by which QFM is shifted
    T_K         float, pd.Series-like
        Temperature in Kelvin


    Returns
    -------
    fO2
    """

    offset = 10 ** logshift

    if not isinstance(Pbar, int):
        T_K = np.array(T_K)
        Pbar = np.array(Pbar)

        QFM_pressure_component = np.zeros(
            shape=[
                len(T_K),
            ]
        )

        for i, (temperature, pressure) in enumerate(zip(T_K, Pbar)):
            QFM_pressure_component[i] = QFM_pressure(temperature, pressure)

    else:
        QFM_pressure_component = QFM_pressure(T_K, Pbar)

    return np.exp((QFM_1bar(T_K) + QFM_pressure_component) / (R * T_K)) * offset


def FeRedox_KC(composition, T_K, fO2, Pbar):
    """
    Calculate Fe-redox equilibrium for silicate liquids according to equation 7 from Kress and Carmichael (1991).

    Parameters
    ----------
    composition     pd.DataFrame
        Liquid major element composition in wt.% oxides
    T_K             float, pd.Series-like
        temperature in Kelvin
    fO2           float, pd.Series-like
        Oxygen fugacity
    P_Pa            float, pd.Series-like
        Pressure in Pascals


    Returns
    -------
    Fe3+/Fe2+ ratio in the liquid

    """
    P_Pa = Pbar * 1e-5

    LNfO2 = np.log(fO2)

    components = ["Al2O3", "FeO", "CaO", "Na2O", "K2O"]

    # Parameters from table 7
    a = 0.196
    b = 1.1492e4
    c = -6.675

    dCoefficients = pd.Series(
        {"Al2O3": -2.243, "FeO": -1.828, "CaO": 3.201, "Na2O": 5.854, "K2O": 6.215},
        name="dCoeff",
    )

    e = -3.36
    f = -7.01e-7
    g = -1.54e-10
    h = 3.85e-17
    T0 = 1673

    molFractions = cc.componentFractions(composition, normalise=True)
    sumComponents = molFractions.loc[:, components].mul(dCoefficients).sum(axis=1)

    part1 = a * LNfO2 + b / T_K + c + sumComponents
    part2 = e * (1 - T0 / T_K - np.log(T_K / T0))
    part3 = f * P_Pa / T_K + g * ((T_K - T0) * P_Pa) / T_K + h * P_Pa ** 2 / T_K

    return np.exp(part1 + part2 + part3)


def FeRedox_Boris(composition: pd.DataFrame, T_K, fO2, *args):
    """
    Borisov et al. (2018), equation 4

    Returns
    -------
    Fe3+/Fe2+ ratio in the liquid
    """   

    oxides = ["SiO2", "TiO2", "MgO", "CaO", "Na2O", "K2O", "Al2O3", "P2O5"]

    missing_oxides = set(oxides).difference(composition.columns)

    if len(missing_oxides) > 0:

        for oxide in missing_oxides:
            composition[oxide] = 0.0

        w.warn(f"{', '.join(str(i) for i in missing_oxides)} missing in composition and set to 0.")

    # Oxide molar fractions
    mol_fractions = cc.componentFractions(composition, normalise=True)

    part1 = (
        0.207 * np.log10(fO2)
        + 4633.3 / T_K
        - 0.445 * mol_fractions["SiO2"]
        - 0.900 * mol_fractions["TiO2"]
        + 1.532 * mol_fractions["MgO"]
    )
    part2 = (
        0.341 * mol_fractions["CaO"]
        + 2.030 * mol_fractions["Na2O"]
        + 3.355 * mol_fractions["K2O"]
        - 4.851 * mol_fractions["P2O5"]
    )
    part3 = (
        -3.081 * mol_fractions["SiO2"] * mol_fractions["Al2O3"]
        - 4.370 * mol_fractions["SiO2"] * mol_fractions["MgO"]
        - 1.852
    )

    return 10**(part1 + part2 + part3)


def FeRedox_QFM(composition, T_K, Pbar, logshift=0, model="Borisov"):
    """
    Calculate Fe-redox equilibrium at QFM oxygen buffer for silicate liquids.
    Uses either equation 7 from Kress and Carmichael (1991) or equation 4 from Borisov et al. (2018).

    Parameters
    ----------
    composition     pd.DataFrame
        Liquid major element composition in wt.% oxides
    T_K             float, pd.Series-like
        temperature in Kelvin
    Pbar            float, pd.Series-like
        Pressure in bars
    logshift        int, pd.Series-like
        log units shift of QFM

    Returns
    -------
    Fe3+/Fe2+ ratio in the liquid

    """

    model_dict = {"KressCarmichael": FeRedox_KC, "Borisov": FeRedox_Boris}
    equation = model_dict[model]

    fO2 = fO2QFM(logshift, T_K, Pbar)

    return equation(composition, T_K, fO2, Pbar)
