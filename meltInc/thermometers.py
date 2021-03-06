from . import chemcalc as cc
import numpy as np
import pandas as pd
from scipy.constants import R #J*K-1*mol-1


def liquid_putirka(liquid, P_bar):

    """Liquid thermometer

    Equation 16 of Putirka (2008) calculates liquiqdus temperature for liquid compositions. Requires equilibrium with olivine + plagioclase + clinopyroxene.

    Parameters
    ----------
    melt : pd.DataFrame
        major element composition of a liquid (wt. %).
    P_bar : int or list-like
        crystallisation pressures. Indices need to match melt if using pd.Series.


    Returns
    -------
    pd.Series
        liquidus temperatures in degrees Kelvin.
    """
    
    if isinstance(P_bar, pd.Series):
        if not liquid.index.equals(P_bar.index):
            raise RuntimeError('Melt and P_bar indices don\'t match')

    oxides = set(['SiO2', 'Al2O3', 'MgO'])
    absentOxides = oxides.difference(liquid.columns)

    if len(absentOxides) > 0:
        raise KeyError(f'{absentOxides} not present in liquid dataframe')

    # Convert pressure from bars to GPa
    P_GPa = P_bar / 1E4

    # Calculate molar oxide fractions
    liquid_mol = cc.componentFractions(liquid, type= 'oxide', normalise= True)

    part_1 = -583 + 3141 * liquid_mol['SiO2'] + 15779 * liquid_mol['Al2O3'] + 1338.6 * liquid_mol['MgO']
    part_2 = - 31440 * liquid_mol['SiO2'] * liquid_mol['Al2O3'] + 77.67 * P_GPa

    return part_1 + part_2 + 273.15


def olivine_putirka(olivine, melt, P_bar, H2O = 0):

    """Olivine liquid thermometer

    Equation 22 of Putirka (2008) calculates crystallisation temperatures for olivine-liquid pairs
    
    Parameters
    ----------
    olivine : pd.DataFrame
        major element compositions of olivines (wt. %).
    melt : pd.DataFrame
        major element compositions of liquids in equilibrium with olivine (wt. %). Indices of olivine and melt dataframes should be equal
    P_bar : int or list-like
        crystallisation pressures for olivines in bar. Indices need to match olivine if using pd.Series.
    H2O : int or list-like
        Water content of melts in the wt.%. Indices need to match olivine if using pd.Series.

    Returns
    -------
    pd.Series
        crystallisation temperatures in degrees Kelvin
    """    

    #convert pressure from bars to GPa
    P_GPa = P_bar / 1E4
    
    if not olivine.index.equals(melt.index):
        raise RuntimeError('olivine and melt indices do not match')
        
    for _,(i,j) in enumerate({'P_bar': P_bar, 'H2O': H2O}.items()):
        if isinstance(j, pd.Series):
            if not olivine.index.equals(j.index):
                raise RuntimeError(f'{i} and olivine indices do not match')     

    #calculate cation fractions
    olivine_cat = cc.componentFractions(olivine, type = 'cation', normalise = True)
    melt_cat = cc.componentFractions(melt, type = 'cation', normalise = True)

    #Mg partitioning
    D_Mg = olivine_cat['MgO']/melt_cat['MgO']

    NM_components = ['FeO', 'MnO', 'MgO', 'CaO', 'CoO', 'NiO']
    C_NM = melt_cat.loc[:, [k for k in NM_components if k in melt_cat.columns]].sum(axis = 1)
    
    NF_components = ['SiO2', 'NaO', 'K2O', 'TiO2']
    C_NF = melt_cat.loc[:,[l for l in NF_components if l in melt_cat.columns]].sum(axis = 1)

    part_1 = 15294.6+1318.8*P_GPa + 2.4834*P_GPa**2.
    part_2 = 8.048 + 2.8352*np.log(D_Mg) + 2.097*np.log(1.5*C_NM)
    part_3 = 1.41*C_NF + 0.222*H2O + 0.5*P_GPa
    
    return (part_1 / (part_2 + part_3)) + 273.15

    
