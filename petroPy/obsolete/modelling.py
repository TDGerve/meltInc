import numpy as np
from scipy.special import erf
from scipy.optimize import fmin


    
def errorfunction(x, amplitude, sign, center, height):
    """Returns an error function evaluated at x"""
        
    return amplitude * erf(sign * (x - center)) + height


    
def errorparams(distance, composition, half_profile= True, uncertainty= 1, xtol= 1E-7, ftol= 1E-7):
    
    if not half_profile:
        center = (max(distance) - min(distance))/2
    elif half_profile:
        center = min(distance)
            
    initial_values = {'a': abs((composition.iloc[-1] - composition[0])/2), #plateau difference composition
                      'b': np.sign(composition.iloc[-1] - composition[0]), # positive/negative slope
                      'c': center, #middle of error function
                      'd': min(composition) + abs((composition.iloc[-1] - composition[0]/2))} #'middle'composition

    
    #Chi squared of profile fited to an error function
    cost_function = lambda x : np.sum(((composition - errorfunction(distance, x[0], x[1], x[2], x[3]))**2)/(uncertainty*2))
    
    optimal_values = fmin(cost_function, x0 = list(initial_values.values()), xtol = xtol, ftol = ftol)
    
    return optimal_values

    
    

