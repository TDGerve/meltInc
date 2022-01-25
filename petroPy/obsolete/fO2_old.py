# ThermoquartzHol_P = [-910710, 41.43]
    # ThermocoesiteHol_P = [-906990, 39.60]
    # ThermostishoviteHol_P = [-876720, 24.00]
    # ThermofayaliteHol_P = [-1477510, 151]  
    # ThermoringwooditeHol_P = [-1471510, 140] 
    # ThermomagnetiteHol_P = [-1114510, 146.9]

    # phases = ['fayalite', 'ringwoodite', 'quartz', 'coesite', 'stishovite', 'magnetite']
    # HollP = {}
    # #volumes could also be calculated in a loop
   
    # #Fayalite
    # Fayaliteth = fO2.teosth(phase= 'fayalite', pkbar= pgpa*10, t= T_K)
    # vfayalite=  10*float(opt.fsolve(lambda v: fO2.teosthV(v, *Fayaliteth[:3], v0= getattr(fO2.EOSparams, 'fayalite')['v0'], pkbar= pgpa*10), 1))

    # ## Ringwoodite
    # Ringwooditeth = fO2.teosth(ThermoringwooditeHol_P[1], v0= 4.203, n= 7, a0= 0.0000222, K= 1977, dKdP= 4.92, dKdP2= -0.0025, pkbar= pgpa*10, tref= 298.15, t= T_K)
    # vringwoodite = 10*float(opt.fsolve(lambda v: fO2.teosthV(v= v, v0= 4.203, a= Ringwooditeth[1], b= Ringwooditeth[2], c= Ringwooditeth[3], pkbar= pgpa*10, pth= Ringwooditeth[0]), 1))

    # ##Quartz_HollP
    # Quartzth = fO2.teosth(s= ThermoquartzHol_P[1], v0= 2.269, n= 3, a0= 0, K= 730, dKdP= 6, dKdP2= -0.0082, pkbar= pgpa*10, tref= 298.15, t= T_K)
    # vquartz = 10*float(opt.fsolve(lambda v: fO2.teosthV(v= v, v0= 2.269, a= Quartzth[1], b= Quartzth[2], c= Quartzth[3], pkbar= pgpa*10, pth= Quartzth[0]), 1))

    # ##Coesite_HollP
    # Coesiteth = fO2.teosth(s= ThermocoesiteHol_P[1], v0= 2.064, n= 3, a0= 0.0000123, K= 979, dKdP= 4.19, dKdP2= -0.0043, pkbar= pgpa*10, tref= 298.15, t= T_K)
    # vcoesite = 10*float(opt.fsolve(lambda v: fO2.teosthV(v= v, v0= 2.064, a= Coesiteth[1], b= Coesiteth[2], c= Coesiteth[3], pkbar= pgpa*10, pth= Coesiteth[0]), 1))

    # ##Stishovite_HollP
    # Stishoviteth = fO2.teosth(s= ThermostishoviteHol_P[1], v0= 1.401, n= 3, a0= 0.0000158, K= 3090, dKdP= 4.6, dKdP2= -0.00150, pkbar= pgpa*10, tref= 298.15, t= T_K)
    # vstishovite = 10*float(opt.fsolve(lambda v: fO2.teosthV(v= v, v0= 1.401, a= Stishoviteth[1], b= Stishoviteth[2], c= Stishoviteth[3], pkbar= pgpa*10, pth= Stishoviteth[0]), 1))

    # ##Magnetite_HollP
    # Magnetiteth = fO2.teosth(ThermomagnetiteHol_P[1], v0= 4.452, n= 7, a0= 3.71e-5, K= 1857, dKdP= 4.05, dKdP2= -0.0022, pkbar= pgpa*10, tref= 298.15, t= T_K)
    # vmagnetite = 10*float(opt.fsolve(lambda v: fO2.teosthV(v= v, v0= 4.452, a= Magnetiteth[1], b= Magnetiteth[2], c= Magnetiteth[3], pkbar= pgpa*10, pth= Magnetiteth[0]), 1))
    