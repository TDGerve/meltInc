from . import chemcalc as cc
import numpy as np
import pandas as pd
from scipy.constants import R  # J*K-1*mol-1
from . import fO2


def Phi_toplis(molar_SiO2, molar_Na2O, molar_K2O):
    """ "returns Phi component for Toplis (2005) Fe-Mg Kd calculations

    Equation 12 from Toplis (2005) calculates a Phi parameter to correct SiO2 for alkali-bearing liquids.
    This expression is only valid for SiO2 <= 60 mol %


    Parameters
    ----------
    molar_SiO2 : int or list-like

    molar_Na2O : int or list-like

    molar_K2O : int or list-like


    Returns
    -------
        int or list-like
    """

    if sum(np.array(molar_SiO2) > 60):
        raise RuntimeError("SiO2 >60 mol% present")

    return (0.46 * (100 / (100 - molar_SiO2)) - 0.93) * (molar_Na2O + molar_K2O) + (
        -5.33 * (100 / (100 - molar_SiO2)) + 9.69
    )


def SiO2_A_toplis(liquid, H2O=None):
    """returns adjusted SiO2 for Toplis (2005) Fe-Mg Kd calculations

    Equations 11 and 14 calculate adjusted molar SiO2 by correcting for akalis and water


    Parameters
    ----------
    molar_SiO2 : int or list-like

    molar_Na2O : int or list-like

    molar_K2O : int or list-like

    Phi : int or list-like
        coefficient for alkali correction, needs to be calculated according to Toplis (eq 12, 2005)

    H2O : int or list-like, optional
        wt. %


    Returns
    -------
    int or list-like
    """

    oxides = cc.oxideweights()
    components = oxides.index.intersection(liquid.columns)

    # Calculate melt molar concentrations
    # Molar fractions normalised to 1
    liq_molar_concentrations = cc.componentFractions(
        liquid, type="oxide", normalise=True
    )
    # Total of the input composition
    total = liquid.loc[:, components].sum(axis=1)
    # Molar fractions renormalised to input total to account for potentially missing volatiles
    liq_molar_concentrations.loc[:, components] = liq_molar_concentrations.loc[
        :, components
    ].mul(total, axis=0)
    liq_molar_concentrations["total"] = liq_molar_concentrations.loc[:, components].sum(
        axis=1
    )

    molar_SiO2 = liq_molar_concentrations["SiO2"]
    molar_Na2O = liq_molar_concentrations["Na2O"]
    molar_K2O = liq_molar_concentrations["K2O"]

    Phi = Phi_toplis(molar_SiO2, molar_Na2O, molar_K2O)
    # Equation 11
    SiO2_A = molar_SiO2 + Phi * (molar_Na2O + molar_K2O)

    if H2O is not None:
        SiO2_A = SiO2_A + 0.8 * H2O  # equation 14

    return SiO2_A


def Kd_toplis(T_K, P_bar, forsterite, SiO2_A):
    # equation 10 of Toplis (2005)

    return np.exp(
        (-6766 / (R * T_K) - 7.34 / R)
        + np.log(0.036 * SiO2_A - 0.22)
        + (3000 * (1 - 2 * forsterite) / (R * T_K))
        + (0.035 * (P_bar - 1) / (R * T_K))
    )


def KdToplis_iterator(
    liquid: pd.DataFrame,
    olivine_forsterite,
    T_K,
    Pbar,
    H2O=None,
    QFMlogshift=0,
    FeRedox_model="Borisov",
    **kwargs
):

    """equation 10 of Toplis (2005) iteratively solved for forsterite content

    Parameters
    ----------
    liquid : pd.DataFrame
        liquid composition in oxide wt.%
    olivine_forsterite
        forsterite fraction in olivine as Mg / (Mg + Fe)
    T_K :
        Temperature in Kelvin
    Pbar :
        Pressure in bar
    H2O : optional
        Water in wt.%
    QFMlogshift : int, default: 0
        log units shifts of fO2 buffered at QFM
    FeRedox_model : str, default: 'Borisov'
        Model to calculate liquid Fe3+/Fe2+, options: 'Borisov', 'KressCarmichael'
    """

    if isinstance(T_K, (int, float)):
        T_K = pd.Series(T_K, index=liquid.index)
    if isinstance(Pbar, (int, float)):
        Pbar = pd.Series(Pbar, index=liquid.index)

    fo_converge_default = 0.001
    fo_converge = kwargs.setdefault("fo_converge", fo_converge_default)

    melts = liquid.copy()
    if (H2O is None) & ("H2O" in melts.columns):
        H2O = melts["H2O"]

    SiO2mol_A = SiO2_A_toplis(liquid, H2O)

    # Fo content observed
    forsterite = olivine_forsterite.copy()
    # Toplis Kd
    melts["Kd"] = Kd_toplis(T_K, Pbar, forsterite, SiO2mol_A)

    # Liquid Fe3/Fe2
    Fe3Fe2 = fO2.FeRedox_QFM(
        liquid, T_K, Pbar, logshift=QFMlogshift, model=FeRedox_model
    )
    # Liquid Fe2+/Fe(total)
    Fe2Fe_total = 1 / (1 + Fe3Fe2)
    # Convert oxides to cations
    oxide_weights = cc.oxideweights()
    ox_fact = oxide_weights["FeO"] / oxide_weights["MgO"]
    # liquid Fe2+/Mg
    Fe2Mg_liquid = (liquid["FeO"] / liquid["MgO"]) / ox_fact * Fe2Fe_total
    # Equilibrium forsterite content according to Kd
    forsterite_EQ = 1 / (1 + melts["Kd"] * Fe2Mg_liquid)

    # Difference between observed Fo and equilibrium Fo
    forsterite_delta = (forsterite - forsterite_EQ) / forsterite

    # iterate until equilibrium forsterite content doesn't change any more
    while sum(forsterite_delta > fo_converge) > 1:

        iterate = forsterite_delta > fo_converge

        melts.loc[iterate, "Kd"] = Kd_toplis(
            T_K[iterate],
            Pbar[iterate],
            forsterite_EQ.loc[iterate],
            SiO2mol_A.loc[iterate],
        )

        forsterite.loc[iterate] = forsterite_EQ.loc[iterate].copy()

        forsterite_EQ.loc[iterate] = 1 / (
            1 + melts.loc[iterate, "Kd"] * Fe2Mg_liquid[iterate]
        )

        forsterite_delta.loc[iterate] = (
            forsterite.loc[iterate] - forsterite_EQ.loc[iterate]
        ) / forsterite.loc[iterate]

    return melts["Kd"]


def Kd_Blundy(forsterite, Fe3Fe2_liquid, T_K):
    """
    Blundy et al., 2020, equation 8
    """

    Fe3FeTotal = Fe3Fe2_liquid / (1 + Fe3Fe2_liquid)

    return 0.3642 * (1 - Fe3FeTotal) * np.exp(312.7 * (1 - 2 * forsterite) / T_K)


def Kd_Blundy_iterator(
    liquid: pd.DataFrame,
    olivine_forsterite,
    T_K,
    Pbar,
    QFMlogshift=0,
    FeRedox_model="Borisov",
    **kwargs
):

    """equation 8 by Blundy (2020) iteratively solved for forsterite content

    Parameters
    ----------
    liquid : pd.DataFrame
        liquid composition in oxide wt.%
    olivine_forsterite
        forsterite fraction in olivine as Mg / (Mg + Fe)
    T_K :
        Temperature in Kelvin
    Pbar :
        Pressure in bar
    QFMlogshift : int, default: 0
        log units shifts of fO2 buffered at QFM
    FeRedox_model : str, default: 'Borisov'
        Model to calculate liquid Fe3+/Fe2+, options: 'Borisov', 'KressCarmichael'
    """

    if isinstance(T_K, (int, float)):
        T_K = pd.Series(T_K, index=liquid.index)
    if isinstance(Pbar, (int, float)):
        Pbar = pd.Series(Pbar, index=liquid.index)

    fo_converge_default = 0.001
    fo_converge = kwargs.setdefault("fo_converge", fo_converge_default)

    melts = liquid.copy()

    # Fo content observed
    forsterite = olivine_forsterite.copy()
    # Liquid Fe3+/Fe2+
    Fe3Fe2_liquid = fO2.FeRedox_QFM(
        melts, T_K, Pbar, logshift=QFMlogshift, model=FeRedox_model
    )
    # Blundy Kd
    melts["Kd"] = Kd_Blundy(forsterite, Fe3Fe2_liquid, T_K)

    # Liquid Fe3/Fe2
    Fe3Fe2 = fO2.FeRedox_QFM(
        liquid, T_K, Pbar, logshift=QFMlogshift, model=FeRedox_model
    )
    # Liquid Fe2+/Fe(total)
    Fe2Fe_total = 1 / (1 + Fe3Fe2)
    # Convert oxides to cations
    oxide_weights = cc.oxideweights()
    ox_fact = oxide_weights["FeO"] / oxide_weights["MgO"]
    # liquid Fe2+/Mg
    Fe2Mg_liquid = (liquid["FeO"] / liquid["MgO"]) / ox_fact * Fe2Fe_total
    # Equilibrium forsterite content according to Kd
    forsterite_EQ = 1 / (1 + melts["Kd"] * Fe2Mg_liquid)

    # Difference between observed Fo and equilibrium Fo
    forsterite_delta = (forsterite - forsterite_EQ) / forsterite

    # iterate until equilibrium forsterite content doesn't change any more
    while sum(forsterite_delta > fo_converge) > 1:

        iterate = forsterite_delta > fo_converge

        melts.loc[iterate, "Kd"] = Kd_Blundy(
            forsterite_EQ.loc[iterate], Fe3Fe2_liquid[iterate], T_K[iterate]
        )

        forsterite.loc[iterate] = forsterite_EQ.loc[iterate].copy()

        forsterite_EQ.loc[iterate] = 1 / (
            1 + melts.loc[iterate, "Kd"] * Fe2Mg_liquid[iterate]
        )

        forsterite_delta.loc[iterate] = (
            forsterite.loc[iterate] - forsterite_EQ.loc[iterate]
        ) / forsterite.loc[iterate]

    return melts["Kd"]


def equilibrium_forsterite(
    liquid, Kd, T_K, Pbar, QFMlogshift=0, FeRedox_model="Borisov"
):
    """
    Parameters
    ----------
    liquid :
        Liquid composition in oxide wt.%
    Kd :
        (Fe_olivine / Fe_liquid) * (Mg_liquid / Mg_olivine) partitioning coefficient
    T_K :
        Temperature in Kelvin
    Pbar :
        Pressure in bar
    QFMlogshift : int, default: 0
        log units shifts of fO2 buffered at QFM
    FeRedox_model : str, default: 'Borisov'
        Model to calculate liquid Fe3+/Fe2+, options: 'Borisov', 'KressCarmichael'

    Returns
    -------
    Equilibrium forsterite fraction as Mg/(Mg + Fe)
    """

    oxide_weights = cc.oxideweights()
    # Convert oxides to cations
    ox_fact = oxide_weights["FeO"] / oxide_weights["MgO"]
    # liquid Fe3/Fe2
    Fe3Fe2 = fO2.FeRedox_QFM(
        liquid, T_K, Pbar, logshift=QFMlogshift, model=FeRedox_model
    )
    # liquid Fe2+/Fe3+ ratio
    Fe2Fe_total = 1 / (1 + Fe3Fe2)
    # liquid Fe2+/Mg
    Fe2Mg_liquid = (liquid["FeO"] / liquid["MgO"]) / ox_fact * Fe2Fe_total

    # Equilibrium forsterite content of olivine according to Kd
    return 1 / (1 + Kd * Fe2Mg_liquid)
