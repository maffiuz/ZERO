# This file contains the builder of the eq. system to give to the solver
from tokamak.vector import vector
import tools.utilities as ut

def Evolution(t: float,values: list) -> list:
    """Main function for the change rate of the various quantities in time

    Args:
        t (float): time of reference
        values (array): the values of the variables

    Returns:
        array: the change rate, evalueted at t (the equations of the manual)
    """
    
    x = vector(values)
    
    # Base DT plasma
    # Ionizations/Recombinations
    rD1 = x.RecIonHandler(ut.REC['D'],Z=1)          # D1 + e- -> D0
    iD0 = x.RecIonHandler(ut.ION['D'],Z=1)          # D0 -> D1 + e-
    rT1 = x.RecIonHandler(ut.REC['T'],Z=1)          # T1 + e- -> T0 
    iT0 = x.RecIonHandler(ut.ION['T'],Z=1)          # T0 -> T1 + e-
    rHe4_2 = x.RecIonHandler(ut.REC['He4'],Z=2)     # He4_2 + e- -> He4_1
    rHe4_1 = x.RecIonHandler(ut.REC['He4'],Z=1)     # He4_1 + e- -> He4_0
    iHe4_1 = x.RecIonHandler(ut.ION['He4'],Z=2)     # He4_1 -> He4_2 + e-
    iHe4_0 = x.RecIonHandler(ut.ION['He4'],Z=1)     # He4_0 -> He4_1 + e-
    
    # Mass losses
    lossD1 = x.MassLoss('D',Z=1)
    lossT1 = x.MassLoss('T',Z=1)
    lossHe4_1 = x.MassLoss('He4',Z=1)
    lossHe4_2 = x.MassLoss('He4',Z=2)
    
    # The lost ions partially recycling back as neutral
    recycling_D0 = x.ParticleRecycling('D',Z=1)
    recycling_T0 = x.ParticleRecycling('T',Z=1)
    recycling_He40 = sum(x.ParticleRecycling('He4',Z) for Z in (1,2))
    
    # Refueling
    refuel_D0 = x.Refueling(t,'D')
    refuel_T0 = x.Refueling(t,'T')

    # Reactions per second
    reactions_second_DT = x.getReactionsSecond('H3(D,N)HE4')
    reactions_second_DD_T = x.getReactionsSecond('H2(D,P)H3')
    if x.nHe3:
        reactions_second_DHE3 = x.getReactionsSecond('HE3(D,P)HE4')
        reactions_second_DD_HE3 = x.getReactionsSecond('H2(D,N)HE3')
    else:
        reactions_second_DHE3 = 0.
        reactions_second_DD_HE3 = 0.
        
    # Building the change rate vector
    change_rate = []
    # Total energy time derivative
    change_rate.append(x.PowerSource(t) - x.PowerSink(t))  
    # D time derivatives
    change_rate.append(rD1-iD0+refuel_D0+recycling_D0)
    change_rate.append(iD0-rD1-reactions_second_DT-2*reactions_second_DD_HE3-2*reactions_second_DD_T-reactions_second_DHE3-lossD1)
    # T time derivatives
    change_rate.append(rT1-iT0+refuel_T0+recycling_T0)
    change_rate.append(iT0-rT1-reactions_second_DT-lossT1+reactions_second_DD_T)
    # He4 time derivatives
    change_rate.append(rHe4_1-iHe4_0+recycling_He40)
    change_rate.append(rHe4_2-iHe4_1+iHe4_0-rHe4_1-lossHe4_1)
    change_rate.append(reactions_second_DT+reactions_second_DHE3-rHe4_2+iHe4_1-lossHe4_2)
    
    # He3 derivatives, activated only in presence of parassitic reactions that involve He3
    if x.nHe3:
        rHe3_2 = x.RecIonHandler(ut.REC['He3'],Z=2)     # He3_2 + e- -> He3_1
        rHe3_1 = x.RecIonHandler(ut.REC['He3'],Z=1)     # He3_1 + e- -> He3_0
        iHe3_1 = x.RecIonHandler(ut.ION['He3'],Z=2)     # He3_1 -> He3_2 + e-
        iHe3_0 = x.RecIonHandler(ut.ION['He3'],Z=1)     # He3_0 -> He3_1 + e-
        
        lossHe3_1 = x.MassLoss('He3',Z=1)
        lossHe3_2 = x.MassLoss('He3',Z=2)
        
        recycling_He30 = sum(x.ParticleRecycling('He3',Z) for Z in (1,2))

        change_rate.append(rHe3_1-iHe3_0+recycling_He30)                                                # He3
        change_rate.append(rHe3_2-iHe3_1+iHe3_0-rHe3_1-lossHe3_1)                                       # He3 +
        change_rate.append(reactions_second_DD_HE3-reactions_second_DHE3-lossHe3_2-rHe3_2+iHe3_1)       # He3 2+
    
    # Impurities change rate
    for element in ut.IMPURITIES_CHARGE_STATES.keys():
        for z in range(ut.IMPURITIES_CHARGE_STATES[element]):
            balance = 0
            # Interaction with lower ionization state
            if z > 0:
                balance -= x.RecIonHandler(ut.REC[element],z)
                balance += x.RecIonHandler(ut.ION[element],z)
                balance -= x.MassLoss('Imp',z,element)
            else:
                balance += x.impurityInjection(t, element)
                balance += sum(x.ParticleRecycling('Imp',Z,element) for Z in range(ut.IMPURITIES_CHARGE_STATES[element]))

            # Interaction with higher ionization state
            if z < ut.IMPURITIES_CHARGE_STATES[element]-1:
                balance -= x.RecIonHandler(ut.ION[element],z+1)
                balance += x.RecIonHandler(ut.REC[element],z+1)

            change_rate.append(balance)
    
    return change_rate      
        