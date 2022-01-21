from . import chemcalc as cc
import numpy as np
import pandas as pd
from scipy.constants import R #J*K-1*mol-1


def Phi_toplis(molar_SiO2, molar_Na2O, molar_K2O):
    """"returns Phi component for Toplis (2005) Fe-Mg Kd calculations

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
        raise RuntimeError('SiO2 >60 mol% present')

    return (0.46*(100/(100-molar_SiO2)) - 0.93) * (molar_Na2O+molar_K2O) + (-5.33*(100/(100-molar_SiO2)) + 9.69)


def SiO2_A_toplis(molar_SiO2, molar_Na2O, molar_K2O, Phi, H2O=None):
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

    SiO2_A = molar_SiO2 + Phi*(molar_Na2O*molar_K2O)  # equation 11

    if H2O:  # is not None
        SiO2_A = SiO2_A + 0.8*H2O  # equation 14

    return SiO2_A


def Kd_toplis(T_K, P_bar, forsterite, SiO2_A):
    # equation 10 of Toplis (2005)

    return np.exp((-6766/(R*T_K) - 7.34/R) + np.log(0.036*SiO2_A - 0.22) + (3000*(1 - 2*forsterite)/(R*T_K)) + (0.035*(P_bar-1)/(R*T_K)))

    
def Kd_ToplisIter(liquid, olivine_Fo, T_K, P, H2O=None):
    # equation 10 of Toplis (2005) iteratively solved for forsterite content
    oxideweights = cc.oxideweights()
    melts = liquid.copy()
    components = [element for element in oxideweights
                  if element in liquid.columns]

    molar_weights = oxideweights.join(pd.DataFrame({'H2O': (1.008*2)+15.999}, index=[0]))
    ox_fact = (molar_weights['MgO']/molar_weights['FeO'])[0]

    # calculate melt into molar fractions. Add water if specified. Normalised to measured totals
    melt_molar = cc.componentFractions(melts)
    if H2O:
        melts = melts.join(H2O)
        melt_molar = melt_molar.join(H2O.divide((1.008*2+15.999)))
        total = melts.loc[:, components].sum(axis=1) + H2O
        components.append('H2O')
    else:
        total = melts.loc[:, components].sum(axis=1)
    melt_molar.loc[:, components] = melt_molar.loc[:, components].mul((total / melt_molar['total']), axis=0)
    melt_molar['total'] = melt_molar.loc[:, components].sum(axis=1)

    # correction of Si according to Toplis (2005)
    Phi = Phi_toplis(melt_molar['SiO2'], melt_molar['Na2O'], melt_molar['K2O'])
    SiO2mol_A = SiO2_A_toplis(melt_molar['SiO2'], melt_molar['Na2O'], melt_molar['K2O'], Phi, H2O)

    forsterite = olivine_Fo.copy()
    melts['Kd'] = Kd_toplis(T_K, P, forsterite, SiO2mol_A)
    forsterite_EQ = 1/(1+ox_fact*melts['Kd']*(melts['FeO']/melts['MgO']))

    #difference between input forsterite content and equilibrium forsterite content according to calculated Kd's
    forsterite_delta = (forsterite - forsterite_EQ)/forsterite

    # iterate until equilibrium forsterite content doesn't change any more
    while sum(forsterite_delta > 0.01) > 1:
        iterate = forsterite_delta > 0.01
        melts.loc[iterate, 'Kd'] = Kd_toplis(T_K, P, forsterite_EQ.loc[iterate], SiO2mol_A.loc[iterate])
        forsterite.loc[iterate] = forsterite_EQ.loc[iterate].copy()
        forsterite_EQ.loc[iterate] = 1/(1+ox_fact*melts.loc[iterate, 'Kd']*(melts.loc[iterate, 'FeO']/melts.loc[iterate, 'MgO']))
        forsterite_delta.loc[iterate] = (forsterite.loc[iterate] - forsterite_EQ.loc[iterate])/forsterite.loc[iterate]

    return melts['Kd']
