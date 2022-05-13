from . import fO2_calc as fO2
from .. import chemcalc as cc
import warnings as w
import numpy as np
import pandas as pd


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

    fO2 = fO2.fO2_QFM(logshift, T_K, Pbar)

    return equation(composition, T_K, fO2, Pbar)
