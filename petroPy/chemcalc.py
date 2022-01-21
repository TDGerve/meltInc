import pandas as pd
import matplotlib.pyplot as plt
from importlib import resources
import numpy as np
import warnings as w
import math
from scipy.constants import R


# %% CHEMISTRY COMPONENTS

def elementweights():
    """Returns dataframe with weights of common elements"""

    return pd.DataFrame({'Si': 28.085,
                         'Al': 26.981,
                         'Mg': 24.305,
                         'Ca': 40.078,
                         'Fe': 55.845,
                         'Na': 22.989,
                         'K': 39.0983,
                         'Mn': 54.938,
                         'Ti': 47.867,
                         'S': 32.06,
                         'P': 30.973,
                         'Cl': 35.45,
                         'Cr': 51.996,
                         'Ni': 58.6934,
                         'Ba': 137.327,
                         'Cu': 63.546,
                         'O': 15.999}, index=[0])


def oxideweights():
    """Returns dataframe with weights of common oxides"""

    element = elementweights()
    O = element['O']

    return pd.DataFrame({'SiO2': (element['Si'] + 2*O),
                         'Al2O3': (2*element['Al'] + 3*O),
                         'MgO': element['Mg'] + O,
                         'CaO': element['Ca'] + O,
                         'FeO': element['Fe'] + O,
                         'Na2O': (2*element['Na'] + O),
                         'K2O': (element['K']*2 + O),
                         'MnO': 54.938 + O,
                         'TiO2': (element['Ti'] + 2*O),
                         'SO3': (element['S'] + 3*O),
                         'P2O5': (element['P']*2 + O*5),
                         'Cl': element['Cl'],
                         'Cr2O3': element['Cr']*2 + 3*O,
                         'NiO': element['Ni']+O,
                         'BaO': element['Ba']+O,
                         'CuO': element['Cu']+O}, index=[0])


def cations():
    """Returns dataframe with number of cations in common oxides"""

    return pd.DataFrame({'SiO2': 1,
                         'Al2O3': 2,
                         'MgO': 1,
                         'CaO': 1,
                         'FeO': 1,
                         'Na2O': 2,
                         'K2O': 2,
                         'MnO': 1,
                         'TiO2': 1,
                         'SO3': 1,
                         'P2O5': 2,
                         'Cl': 1,
                         'Cr2O3': 2,
                         'NiO': 1,
                         'BaO': 1,
                         'CuO': 1}, index=[0])


def oxygens():
    """Returns a dataframe with number of oxygens in common oxides"""

    return pd.DataFrame({'SiO2': 2,
                         'Al2O3': 1.5,
                         'MgO': 1,
                         'CaO': 1,
                         'FeO': 1,
                         'Na2O': 0.5,
                         'K2O': 0.5,
                         'MnO': 1,
                         'TiO2': 2,
                         'SO3': 3,
                         'P2O5': 2.5,
                         'Cl': 0,
                         'Cr2O3': 3,
                         'NiO': 1,
                         'BaO': 1,
                         'CuO': 1}, index=[0])


def componentFractions(composition, type='oxide', normalise=False, normFactor=None, elements=None):
    """Calulate oxide or cation fractions from major element compositions, with optional 
    normalisation to a fixed number of oxygens,cations or the cation,oxide total

    Parameters
    ----------
    composition : pandas.DataFrame
        pandas dataframe with major element composition in oxide wt.%
    type : {'oxide', 'cation'}, default: 'oxide
        component type
    normalise : {False, 'O', 'cat', 'total'}, optional
        normalisation type, False for no normalisation, 'O' for oxygen
        'cat' for cation and 'total' for total amount of components (oxide 
        or cation - specified in 'type')
    normFactor : int, optional
        amount of oxides or cations to normalise to, e.g. 8 oxygen for plagioclase
    elements : list, optional
        elements to use in the calculations

    Returns
    -------
    pandas.DataFrame
        pandas dataframe with component fractions, cation totals and oxygen totals
    """

    if type not in ['oxide', 'cation']:
        raise ValueError('type should be \'oxide\' or \'cation\'')
    if normalise not in [False, 'O', 'cat', 'total']:
        raise ValueError('type should be False, \'O\', \'cat\' or \'total\'')

    if isinstance(composition, pd.Series):
        composition = pd.DataFrame(composition).T

    if elements is not None:
        components = elements.copy()
    else:
        components = list(oxideweights().keys())
    components = [x for x in components if x in composition.columns]

    compositionMol = composition.loc[:, components].div(oxideweights().loc[0, components])

    if type == 'cation':
        compositionMol = compositionMol.mul(cations().loc[0, components])

        if normalise == 'O':
            oxygen = compositionMol.loc[:, components].mul(oxygens().loc[0, components], axis=1)
            oxygen['total'] = oxygen.sum(axis=1)

            compositionMol.loc[:, components] = compositionMol.loc[:, components].div(oxygen['total'], axis=0) * normFactor

            compositionMol['cations'] = compositionMol.loc[:,components].sum(axis=1)
            compositionMol['O'] = compositionMol.loc[:, components].mul(oxygens().loc[0, components], axis=1).sum(axis=1)

        if normalise == 'cat':
            compositionMol['cations'] = compositionMol.loc[:,components].sum(axis=1)

            compositionMol.loc[:, components] = compositionMol.loc[:, components].div(compositionMol['cations'], axis=0) * normFactor

            compositionMol['cations'] = compositionMol.loc[:,components].sum(axis=1)
            compositionMol['O'] = compositionMol.loc[:, components].mul(oxygens().loc[0, components], axis=1).sum(axis=1)

        catDict = {i: j for i, j in zip(list(oxideweights().keys()), list(elementweights().keys())[:-1])}
        compositionMol.rename(columns=catDict, inplace=True)

    if normalise == 'total':
        compositionMol['total'] = compositionMol.loc[:, components].sum(axis=1)
        compositionMol.loc[:, components] = compositionMol.loc[:, components].div(compositionMol['total'], axis=0)
        compositionMol['total'] = compositionMol.loc[:, components].sum(axis=1)

    return compositionMol


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

    pxFormula = componentFractions(
        composition, type='cation', normalise='O', normFactor=6)

    Jd = np.zeros(pxFormula.shape[0])
    CaTs = np.zeros(pxFormula.shape[0])
    CaTi = np.zeros(pxFormula.shape[0])
    CrCaTs = np.zeros(pxFormula.shape[0])
    DiHd = np.zeros(pxFormula.shape[0])
    EnFs = np.zeros(pxFormula.shape[0])

    Al_IV = 2 - pxFormula['Si']
    Al_VI = pxFormula['Al'] - Al_IV

    # Jadeite
    Jd[Al_VI < pxFormula['Na']] = Al_VI[Al_VI < pxFormula['Na']]
    Jd[pxFormula['Na'] < Al_VI] = pxFormula['Na'][pxFormula['Na'] < Al_VI]

    # Ca-Tschermak
    CaTs[(Al_VI - Jd) > 0] = Al_VI[(Al_VI - Jd) > 0] - Jd[(Al_VI - Jd) > 0]

    # Ca-Ti
    CaTi[Al_IV > CaTs] = (Al_IV[Al_IV > CaTs] - CaTs[Al_IV > CaTs]) / 2

    # Cr-Ca-Tschermak
    CrCaTs = pxFormula['Cr'].values / 2

    # Diopside-Hedenbergite
    DiHd = pxFormula['Ca'].values - CaTi - CaTs - CrCaTs

    EnFs = (pxFormula['Fe'].values + pxFormula['Mg'].values - DiHd) / 2

    components = pd.DataFrame({'Jd': Jd,
                              'CaTs': CaTs,
                               'CaTi': CaTi,
                               'CrCaTs': CrCaTs,
                               'DiHd': DiHd,
                               'EnFs': EnFs
                               })

    return components


# %% RESERVOIR COMPOSITIONS

def C1chondrite():
    """Returns C1 chondrite composition from McDonough & Sun (1995)"""

    with resources.open_text('petroPy.static', 'Mcdonough_sun_1995.csv') as df:
        C1raw = pd.read_csv(df)

    C1 = C1raw[['C1']].transpose()
    C1.columns = C1raw['Element']

    return C1


def primitiveMantle():
    """Returns primitive mantle composition from McDonough & Sun (1995)"""

    with resources.open_text('petroPy.static', 'Mcdonough_sun_1995.csv') as df:
        PMraw = pd.read_csv(df)

    PM = PMraw[['Pyrolite']].transpose()
    PM.columns = PMraw['Element']

    return PM

# %%


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

    if valency not in [2, 3, 'REE']:
        raise Exception("valency should be 2, 3 or 'REE'")
        return

    valency = {2: 0, 3: 1, 'REE': 1}[valency]

    divalent = {'Mg': 0.89, 'Ba': 1.42, 'Ca': 1.12, 'Eu': 1.25, 'Sr': 1.26}

    REE = {'La': 1.16,
           'Ce': 1.143,
           'Pr': 1.126,
           'Nd': 1.109,
           'Sm': 1.079,
           'Eu': 1.066,
           'Gd': 1.053,
           'Tb': 1.040,
           'Dy': 1.027,
           'Y': 1.019,
           'Ho': 1.015,
           'Er': 1.004,
           'Tm': 0.994,
           'Yb': 0.985,
           'Lu': 0.977
           }

    total = [divalent, REE]

    return total[valency]

# %%


def TAS():
    """Returns a line plot element of classification of volcanic rocks
    in total-alkali vs silica plots
    """

    with resources.open_text('petroPy.static', 'TAS.csv') as df:
        TAS = pd.read_csv(df)

    # TAS= pd.read_csv('D:/Dropbox/python/packages/petroPy/TAS.csv')

    for id in TAS.id.unique():
        plt.plot(TAS.loc[TAS.id == id, 'x'],
                 TAS.loc[TAS.id == id, 'y'], color='k')



