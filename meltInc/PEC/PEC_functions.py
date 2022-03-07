from .. import Kd as kd
from .. import fO2
from .. import chemcalc as cc
import pandas as pd
import numpy as np
from typing import Sequence


def FeMg_olivine_liquid_Kd(
    liquid: pd.DataFrame,
    olivine: pd.DataFrame,
    T_K,
    Pbar,
    model="Blundy",
    QFMlogshift=0,
    FeRedox_model="Borisov",
    **kwargs
):
    """
    Parameters
    ----------
    liquid : pd.DataFrame
        liquid composition in oxide wt.%
    olivine : pd.DataFrame
        olivine  composition in oxide wt.%
    T_K :
        Temperature in Kelvin
    Pbar :
        Temperature in bar
    model : str
        Kd model selection, options: 'Toplis', 'Blundy'

    Other Parameters
    ----------------
    **kwargs :  QFMlofshift, default: 0
                FeRedox_model, default: 'Borisov'

    Returns
    -------

    """

    H2Odefault = None
    H2O = kwargs.setdefault("H2O", H2Odefault)

    # Oxide wt.% to cation fractions
    liquid_cation_fractions = cc.componentFractions(liquid, "cation", normalise=True)
    olivine_cation_fractions = cc.componentFractions(
        olivine, "cation", normalise=True
    )

    # liquid Fe3/Fe2
    Fe3Fe2 = fO2.FeRedox_QFM(
        liquid, T_K, Pbar, logshift=QFMlogshift, model=FeRedox_model
    )

    # liquid Fe2+/Fe3+ ratio
    Fe2Fe_total = 1 / (1 + Fe3Fe2)

    Kd_observed = (
        (olivine_cation_fractions["Fe"] / (liquid_cation_fractions["Fe"] * Fe2Fe_total))
        * (liquid_cation_fractions["Mg"] / olivine_cation_fractions["Mg"])
    ).rename("Kd_observed")

    # Forsterite content of olivines
    Fo = olivine_cation_fractions["Mg"] / (
        olivine_cation_fractions["Fe"] + olivine_cation_fractions["Mg"]
    )

    if model == "Toplis":
        Kd_calculated = kd.KdToplis_iterator(
            liquid,
            Fo,
            T_K,
            Pbar,
            H2O=H2O,
            QFMlogshift=QFMlogshift,
            FeRedox_model=FeRedox_model,
        )
    elif model == "Blundy":
        Kd_calculated = kd.Kd_Blundy(Fo, Fe3Fe2, T_K)
    Kd_calculated.rename("Kd_calculated", inplace=True)

    return Kd_observed, Kd_calculated


def FeMg_OlLiq_disequilibrium(Kd_observed, Kd_calculated, **kwargs):

    toleranceDefault = 0.01
    Kd_tolerance = kwargs.setdefault("Kd_tolerance", toleranceDefault)

    return pd.Series(
        np.invert(np.isclose(np.array(Kd_observed), Kd_calculated, atol=Kd_tolerance)),
        index=Kd_observed.index,
    )
