import pandas as pd
import matplotlib.pyplot as plt
from importlib import resources
import numpy as np
import warnings as w
import math
from scipy.constants import R
from typing import List


def elementweights():
    """Returns a pd.Series with weights of common elements"""

    return pd.Series(
        {
            "Si": 28.085,
            "Al": 26.981,
            "Mg": 24.305,
            "Ca": 40.078,
            "Fe": 55.845,
            "Na": 22.989,
            "K": 39.0983,
            "Mn": 54.938,
            "Ti": 47.867,
            "S": 32.06,
            "P": 30.973,
            "Cl": 35.45,
            "Cr": 51.996,
            "Ni": 58.6934,
            "Ba": 137.327,
            "Cu": 63.546,
            "H": 1.008,
            "O": 15.999,
        }
    )


def oxideweights():
    """Returns a pd.Series with weights of common oxides"""

    element = elementweights()
    O = element["O"]

    return pd.Series(
        {
            "SiO2": element["Si"] + 2 * O,
            "Al2O3": 2 * element["Al"] + 3 * O,
            "MgO": element["Mg"] + O,
            "CaO": element["Ca"] + O,
            "FeO": element["Fe"] + O,
            "Na2O": 2 * element["Na"] + O,
            "K2O": element["K"] * 2 + O,
            "MnO": element["Mn"] + O,
            "TiO2": element["Ti"] + 2 * O,
            "SO3": element["S"] + 3 * O,
            "P2O5": element["P"] * 2 + O * 5,
            "Cl": element["Cl"],
            "Cr2O3": element["Cr"] * 2 + 3 * O,
            "NiO": element["Ni"] + O,
            "BaO": element["Ba"] + O,
            "CuO": element["Cu"] + O,
            "H2O": element["H"] * 2 + O,
        }
    )


def cations():
    """Returns a pd.Series with number of cations in common oxides"""

    return pd.Series(
        {
            "SiO2": 1,
            "Al2O3": 2,
            "MgO": 1,
            "CaO": 1,
            "FeO": 1,
            "Na2O": 2,
            "K2O": 2,
            "MnO": 1,
            "TiO2": 1,
            "SO3": 1,
            "P2O5": 2,
            "Cl": 1,
            "Cr2O3": 2,
            "NiO": 1,
            "BaO": 1,
            "CuO": 1,
            "H2O": 1
        }
    )


def oxygens():
    """Returns a pd.Series with number of oxygen per cation in common oxides"""

    return pd.Series(
        {
            "SiO2": 2,
            "Al2O3": 1.5,
            "MgO": 1,
            "CaO": 1,
            "FeO": 1,
            "Na2O": 0.5,
            "K2O": 0.5,
            "MnO": 1,
            "TiO2": 2,
            "SO3": 3,
            "P2O5": 2.5,
            "Cl": 0,
            "Cr2O3": 1.5,
            "NiO": 1,
            "BaO": 1,
            "CuO": 1,
            "H2O": 2
        }
    )


def componentFractions(
    composition: pd.DataFrame,
    type: str = "oxide",
    normalise = False,
    normFactor: int = 1,
    elements: List[str] = None,
):

    """Calulate oxide or cation fractions from major element compositions, with optional
    normalisation to a fixed number of oxygens,cations or the cation,oxide total

    Parameters
    ----------
    composition : pandas.DataFrame
        pandas dataframe with major element composition in oxide wt.%
    type : {'oxide', 'cation'}, default: 'oxide
        component type
    normalise : {False, 'total', 'cation', 'O'}, optional
        normalisation type, False for no normalisation, True for oxide of cation (as set in 'type') and 'O' for oxygen.
    normFactor : int, optional
        amount of oxides, cations or oxygen to normalise to, e.g. 8 oxygen for plagioclase, 100 for concentrations and 1 for fractions
    elements : list, optional
        elements to use in the calculations

    Returns
    -------
    pandas.DataFrame
        pandas dataframe with component fractions, cation totals and oxygen totals
    """

    if type not in ["oxide", "cation"]:
        raise ValueError("type should be 'oxide' or 'cation'")
    if normalise not in [True, False, "total", "cation", "O"]:
        raise ValueError("normalise should be True, False, 'total', 'cation', or 'O'")

    if isinstance(composition, pd.Series):
        composition = pd.DataFrame(composition).T

    if elements is not None:
        components = elements.copy()
    else:
        components = list(oxideweights().index)
    components = composition.columns.intersection(components)

    molar_proportions = composition.loc[:, components].div(oxideweights()[components])

    if type == "oxide":

        if normalise in [True, "total"]:
            molar_proportions["total"] = molar_proportions.loc[:, components].sum(
                axis=1
            )
            molar_proportions.loc[:, components] = (
                molar_proportions.loc[:, components].div(
                    molar_proportions["total"], axis=0
                )
                * normFactor
            )

        molar_proportions["total"] = molar_proportions.loc[:, components].sum(axis=1)

        return molar_proportions

    if type == "cation":

        cation_proportions = molar_proportions.mul(cations()[components])

        oxygen = cation_proportions.loc[:, components].mul(oxygens()[components])
        cation_proportions["O"] = oxygen.sum(axis=1)

        if normalise == "total":

            components_O = np.append(components, "O")

            cation_proportions["total"] = cation_proportions.loc[:, components_O].sum(
                axis=1
            )
            cation_proportions.loc[:, components_O] = (
                cation_proportions.loc[:, components_O].div(
                    cation_proportions["total"], axis=0
                )
                * normFactor
            )
            cation_proportions["total"] = cation_proportions.loc[:, components_O].sum(
                axis=1
            )

        if normalise == "O":

            cation_proportions.loc[:, components] = (
                cation_proportions.loc[:, components].div(
                    cation_proportions["O"], axis=0
                )
                * normFactor
            )

            cation_proportions["cations"] = cation_proportions.loc[:, components].sum(
                axis=1
            )
            cation_proportions["O"] = (
                cation_proportions.loc[:, components]
                .mul(oxygens()[components])
                .sum(axis=1)
            )

        if normalise == "cation":
            cation_proportions["cations"] = cation_proportions.loc[:, components].sum(
                axis=1
            )

            cation_proportions.loc[:, components] = (
                cation_proportions.loc[:, components].div(
                    cation_proportions["cations"], axis=0
                )
                * normFactor
            )

            cation_proportions["cations"] = cation_proportions.loc[:, components].sum(
                axis=1
            )
            cation_proportions["O"] = (
                cation_proportions.loc[:, components]
                .mul(oxygens()[components], axis=1)
                .sum(axis=1)
            )

        catDict = {
            i: j
            for i, j in zip(
                list(oxideweights().index), list(elementweights().index)[:-1]
            )
        }
        cation_proportions.rename(columns=catDict, inplace=True)

        return cation_proportions


def pyroxeneComponents(composition):
    """Calculate components from pyroxene major element compositions

    Returns jadeite, Ca-Tschermak, Ca-Ti, Cr-Ca-Tschermak, Diopside-Hedenbergite and
    Enstatite-Ferrosillite components for pyroxene major element compositions.

    Parameters
    ----------
    composition : pandas.DataFrame
        Pandas dataframe with pyroxene major element composition in oxide wt.%

    Returns
    -------
    pandas.DataFrame
        Pandas dataframe with pyroxene components
    """

    # Calculate cation fractions normalised to 6 oxygen
    pxFormula = componentFractions(
        composition, type="cation", normalise="O", normFactor=6
    )

    Jd = np.zeros(pxFormula.shape[0])
    CaTs = np.zeros(pxFormula.shape[0])
    CaTi = np.zeros(pxFormula.shape[0])
    CrCaTs = np.zeros(pxFormula.shape[0])
    DiHd = np.zeros(pxFormula.shape[0])
    EnFs = np.zeros(pxFormula.shape[0])

    Al_IV = 2 - pxFormula["Si"]
    Al_VI = pxFormula["Al"] - Al_IV

    # Jadeite
    Jd[Al_VI < pxFormula["Na"]] = Al_VI[Al_VI < pxFormula["Na"]]
    Jd[pxFormula["Na"] < Al_VI] = pxFormula["Na"][pxFormula["Na"] < Al_VI]

    # Ca-Tschermak
    CaTs[(Al_VI - Jd) > 0] = Al_VI[(Al_VI - Jd) > 0] - Jd[(Al_VI - Jd) > 0]

    # Ca-Ti
    CaTi[Al_IV > CaTs] = (Al_IV[Al_IV > CaTs] - CaTs[Al_IV > CaTs]) / 2

    # Cr-Ca-Tschermak
    CrCaTs = pxFormula["Cr"].values / 2

    # Diopside-Hedenbergite
    DiHd = pxFormula["Ca"].values - CaTi - CaTs - CrCaTs

    EnFs = (pxFormula["Fe"].values + pxFormula["Mg"].values - DiHd) / 2

    components = pd.DataFrame(
        {
            "Jd": Jd,
            "CaTs": CaTs,
            "CaTi": CaTi,
            "CrCaTs": CrCaTs,
            "DiHd": DiHd,
            "EnFs": EnFs,
        }
    )

    return components


def C1chondrite():
    """Returns C1 chondrite composition from McDonough & Sun (1995)"""

    with resources.open_text("meltInc.static", "Mcdonough_sun_1995.csv") as df:
        C1raw = pd.read_csv(df)

    C1 = C1raw[["C1"]].transpose()
    C1.columns = C1raw["Element"]

    return C1


def primitiveMantle():
    """Returns primitive mantle composition from McDonough & Sun (1995)"""

    with resources.open_text("meltInc.static", "Mcdonough_sun_1995.csv") as df:
        PMraw = pd.read_csv(df)

    PM = PMraw["Pyrolite"].transpose()
    PM.index = PMraw["Element"]

    return PM


def radii(valency):
    """Returns a dictionary with radii for common trace elements

    Parameters
    ----------
    valency : {2, 3, 'REE'}
        valency of requisted trace elements

    Returns
    -------
    dictionary
        dictionary with ionic radii in Angstrom
    """

    if valency not in [2, 3, "REE"]:
        raise Exception("valency should be 2, 3 or 'REE'")
        return

    valency = {2: 0, 3: 1, "REE": 1}[valency]

    divalent = pd.Series({"Mg": 0.89, "Ba": 1.42, "Ca": 1.12, "Eu": 1.25, "Sr": 1.26})

    REE = pd.Series(
        {
            "La": 1.16,
            "Ce": 1.143,
            "Pr": 1.126,
            "Nd": 1.109,
            "Sm": 1.079,
            "Eu": 1.066,
            "Gd": 1.053,
            "Tb": 1.040,
            "Dy": 1.027,
            "Y": 1.019,
            "Ho": 1.015,
            "Er": 1.004,
            "Tm": 0.994,
            "Yb": 0.985,
            "Lu": 0.977,
        }
    )

    total = [divalent, REE]

    return total[valency]
