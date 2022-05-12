from . import eos
from .. import chemcalc as cc
import scipy.optimize as opt
import warnings as w
import numpy as np
import pandas as pd
from scipy.constants import R
import itertools as it


def QFM_pressure_phaseTransitions(T_K, Pbar):
    """
    Pressure contribution to the chemical potential of oxygen (fO2) at QFM
    """

    p_kbar = Pbar / 1e3

    # Pressure of transition SiO2 polymorphs
    qtz_coe = lambda P: eos.phaseTransition(
        P, t=T_K, phase_1="quartz", phase_2="coesite"
    )
    P_qtz_coe = opt.fsolve(qtz_coe, 8)

    coe_stish = lambda P: eos.phaseTransition(
        P, t=T_K, phase_1="coesite", phase_2="stishovite"
    )
    P_coe_stish = opt.fsolve(coe_stish, 8)

    # Pressure of transition Fe2SiO4 polymorphs
    fay_ring = lambda P: eos.phaseTransition(
        P, t=T_K, phase_1="fayalite", phase_2="ringwoodite"
    )
    P_fay_ring = opt.fsolve(fay_ring, 8)

    # VdP SiO2 polymorphs
    *_, VdP_qtz = eos.tait_eos_pressure(phase="quartz", pkbar=min(p_kbar, P_qtz_coe), t=T_K)
    Gibbs_landau = eos.landau_P_dependent(phase="quartz", pkbar=min(p_kbar, P_qtz_coe), t=T_K)
    VdP_qtz = VdP_qtz + Gibbs_landau
    if p_kbar > P_qtz_coe:
        VdP_coe = (
            eos.tait_eos_pressure(phase="coesite", pkbar=min(p_kbar, P_coe_stish), t=T_K)[4]
            - eos.tait_eos_pressure(phase="coesite", pkbar=P_qtz_coe, t=T_K)[4]
        )
        VdP_qtz = VdP_qtz + VdP_coe
        if p_kbar > P_coe_stish:
            VdP_stish = (
                eos.tait_eos_pressure(phase="stishovite", pkbar=p_kbar, t=T_K)[4]
                - eos.tait_eos_pressure(phase="stishovite", pkbar=P_coe_stish, t=T_K)[4]
            )
            VdP_qtz = VdP_qtz + VdP_stish

    # VdP Fe2SiO4 polymorphs
    *_, VdP_fay = eos.tait_eos_pressure(phase="fayalite", pkbar=min(p_kbar, P_fay_ring), t=T_K)
    VdP_Fe2SiO4 = VdP_fay
    if p_kbar > P_fay_ring:
        VdP_ring = (
            eos.tait_eos_pressure(phase="ringwoodite", pkbar=p_kbar, t=T_K)[4]
            - eos.tait_eos_pressure(phase="ringwoodite", pkbar=P_fay_ring, t=T_K)[4]
        )
        VdP_Fe2SiO4 = VdP_Fe2SiO4 + VdP_ring

    # VdP magnetite
    VdP_mt = eos.tait_eos_pressure(phase="magnetite", pkbar=p_kbar, t=T_K)[4]

    return 3e3 * VdP_qtz + 2e3 * VdP_mt - 3e3 * VdP_Fe2SiO4




def QFM_VdP(T_K, Pbar):

    p_kbar = Pbar / 1e3

    # VdP SiO2
    *_, VdP_qtz = eos.tait_eos_pressure(phase="quartz", pkbar=p_kbar, t=T_K)
    Gibbs_landau = eos.landau_P_dependent(phase="quartz", pkbar=p_kbar, t=T_K)
    VdP_qtz = VdP_qtz + Gibbs_landau

    *_, VdP_fay = eos.tait_eos_pressure(phase="fayalite", pkbar=p_kbar, t=T_K)

    *_, VdP_mt = eos.tait_eos_pressure(phase="magnetite", pkbar=p_kbar, t=T_K)

    return VdP_qtz, VdP_mt, VdP_fay


def QFM_VdP_phaseTransitions(T_K, Pbar):

    VdP_SiO2, VdP_mt, VdP_Fe2SiO4 = QFM_VdP(T_K, Pbar)

    p_kbar = Pbar / 1e3

        # Pressure of transition SiO2 polymorphs
    qtz_coe = lambda P: eos.phaseTransition(
        P, t=T_K, phase_1="quartz", phase_2="coesite"
    )
    P_qtz_coe = opt.fsolve(qtz_coe, 8)

    coe_stish = lambda P: eos.phaseTransition(
        P, t=T_K, phase_1="coesite", phase_2="stishovite"
    )
    P_coe_stish = opt.fsolve(coe_stish, 8)

    # Pressure of transition Fe2SiO4 polymorphs
    fay_ring = lambda P: eos.phaseTransition(
        P, t=T_K, phase_1="fayalite", phase_2="ringwoodite"
    )
    P_fay_ring = opt.fsolve(fay_ring, 8)

    if p_kbar > P_qtz_coe:
        VdP_coe = (
            eos.tait_eos_pressure(phase="coesite", pkbar=min(p_kbar, P_coe_stish), t=T_K)[4]
            - eos.tait_eos_pressure(phase="coesite", pkbar=P_qtz_coe, t=T_K)[4]
        )
        VdP_SiO2 = VdP_SiO2 + VdP_coe
        if p_kbar > P_coe_stish:
            VdP_stish = (
                eos.tait_eos_pressure(phase="stishovite", pkbar=p_kbar, t=T_K)[4]
                - eos.tait_eos_pressure(phase="stishovite", pkbar=P_coe_stish, t=T_K)[4]
            )
            VdP_SiO2 = VdP_SiO2 + VdP_stish

    # VdP Fe2SiO4 polymorphs
    if p_kbar > P_fay_ring:
        VdP_ring = (
            eos.tait_eos_pressure(phase="ringwoodite", pkbar=p_kbar, t=T_K)[4]
            - eos.tait_eos_pressure(phase="ringwoodite", pkbar=P_fay_ring, t=T_K)[4]
        )
        VdP_Fe2SiO4 = VdP_Fe2SiO4 + VdP_ring
    
    return VdP_SiO2, VdP_mt, VdP_Fe2SiO4


def QFM_muO2_P(T_K, Pbar):
    """
    Calculate the chemical potential of O2 at pressure P and temperature T from the reaction:

    3 * Fe2SiO4 + O2 = 3 * SiO2 + 2 * Fe3O4

    Parameters
    ----------
    T_K     int, float
        Temperature in Kelvin
    Pbar    int, float
        Pressure in bar

    Returns
    -------
    float
        Chemical potential of O2    
    """


    VdP_quartz, VdP_magnetite, VdP_fayalite = QFM_VdP(T_K, Pbar)
    #kiloJoule to Joule
    muO2 = 1e3 * (3 * VdP_quartz + 2 * VdP_magnetite - 3 * VdP_fayalite)
    
    return muO2


def QFM_muO2_P_phaseTransitions(T_K, Pbar):

    Pbar_is_int = isinstance(Pbar, (int, float))
    T_K_is_int = isinstance(T_K, (int, float))

    # If P and T are not both single numbers
    if not (Pbar_is_int and T_K_is_int):
    
        # If only one variable, P or T, is integer-like
        if bool(Pbar_is_int) ^ bool(T_K_is_int):

            # Cycle the short variable
            T_K = [np.array(T_K), it.cycle(np.array([T_K]))][T_K_is_int]
            Pbar = [np.array(Pbar), it.cycle(np.array([Pbar]))][Pbar_is_int]


        muO2 = np.zeros(shape=[len(T_K),])

        for i, (temperature, pressure) in enumerate(zip(T_K, Pbar)):
            VdP_quartz, VdP_magnetite, VdP_fayalite = QFM_VdP_phaseTransitions(temperature, pressure)
            #kiloJoule to Joule
            muO2[i] = 1e3 * (3 * VdP_quartz + 2 * VdP_magnetite - 3 * VdP_fayalite)

    else:
        VdP_quartz, VdP_magnetite, VdP_fayalite = QFM_VdP_phaseTransitions(T_K, Pbar)
        muO2 = 1e3 * (3 * VdP_quartz + 2 * VdP_magnetite - 3 * VdP_fayalite)

    return muO2  



def QFM_1bar(T_K):
    """
    calculate chemical potential of oxygen at QFM a 1 bar. Equation from O'Neill 1987

    Parameters
    ----------
    T_K     list-like, float
        Temperature in Kelvin

    Returns
    -------
    Chemical potential
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


def fO2QFM_1bar(T_K, logshift=0):
    """
    calculate fO2 at QFM + logshift a 1 bar. Equation from O'Neill 1987

    Parameters
    ----------
    T_K     list-like, float
        Temperature in Kelvin
    logshift   int, float
        Log units by which QFM is shifted

    Returns
    -------
    fO2
    """
    mu_O2 = QFM_1bar(T_K)

    offset = 10 ** logshift

    return np.exp(mu_O2 / (R * T_K)) * offset

def fO2QFM_new(logshift, T_K, Pbar, phaseTransitions=False):
    """
    
    """
    offset = 10 ** logshift

    if not phaseTransitions:
        muO2_pressure = QFM_muO2_P(T_K, Pbar)
    else:
        muO2_pressure = QFM_muO2_P_phaseTransitions(T_K, Pbar)

    muO2_1bar = QFM_1bar(T_K)

    return np.exp((muO2_1bar + muO2_pressure) / (R * T_K)) * offset


def fO2QFM(logshift, T_K, Pbar):
    """
    calculate oxygen fugacity at QFM offset by arbitraty log units.

    Parameters
    ----------
    logshift   int, float
        Log units by which QFM is shifted
    T_K         float, pd.Series-like
        Temperature in Kelvin


    Returns
    -------
    fO2
    """

    offset = 10 ** logshift

    Pbar_is_int = isinstance(Pbar, (int, float))
    T_K_is_int = isinstance(T_K, (int, float))

    # If P and T are not both integer-like
    if not (Pbar_is_int and T_K_is_int):
    
        # If only one variable, P or T, is integer-like
        if bool(Pbar_is_int) ^ bool(T_K_is_int):

            # Cycle the short variable
            T_K = [np.array(T_K), it.cycle(np.array([T_K]))][T_K_is_int]
            Pbar = [np.array(Pbar), it.cycle(np.array([Pbar]))][Pbar_is_int]


        QFM_pressure_component = np.zeros(
            shape=[
                len(T_K),
            ]
        )

        for i, (temperature, pressure) in enumerate(zip(T_K, Pbar)):
            QFM_pressure_component[i] = QFM_pressure_phaseTransitions(temperature, pressure)

    else:
        QFM_pressure_component = QFM_pressure_phaseTransitions(T_K, Pbar)

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
    part3 = f * P_Pa / T_K + g * ((T_K - T0) * P_Pa) / T_K + h * P_Pa**2 / T_K

    return 2 * np.exp(part1 + part2 + part3)


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
    composition pd.DataFrame
        Liquid major element composition in wt.% oxides
    T_K         float, pd.Series-like
        temperature in Kelvin
    Pbar        float, pd.Series-like
        Pressure in bars
    logshift    int, pd.Series-like
        log units shift of QFM
    model       string
        'KressCarmichael' or 'Borisov'       

    Returns
    -------
    Fe3+/Fe2+ ratio in the liquid

    """

    model_dict = {"KressCarmichael": FeRedox_KC, "Borisov": FeRedox_Boris}
    equation = model_dict[model]

    fO2 = fO2QFM(logshift, T_K, Pbar)

    return equation(composition, T_K, fO2, Pbar)
