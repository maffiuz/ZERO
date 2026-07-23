import matplotlib.pyplot as plt
import numpy as np
import userDefined.user_inputs as inp
import tools.utilities as ut
from tokamak.vector import vector

def plotter():    
    # Fetching the data
    base_data = np.loadtxt("out/simulator_data",skiprows=3)

    for element in inp.IMPURITIES.keys():
        path = f'out/impurities_{element}_data'
        try:
            data = np.loadtxt(path,skiprows=5)
            base_data = np.hstack((base_data, data[:,1:]))
        except:
            raise Exception(f'Data for {element} impurity not found')

    # Definition of the plotted quantities
    Density_D = [base_data[:,2],base_data[:,3]]
    Density_He3 = [base_data[:,9],base_data[:,10],base_data[:,11]] if inp.REACTIONS['H2(D,N)HE3'] or inp.REACTIONS['HE3(D,P)HE4'] else []
    Density_He4 = [base_data[:,6],base_data[:,7],base_data[:,8]]
    Density_Impurities = [[] for _ in range(len(inp.IMPURITIES.keys()))]
    Density_T = [base_data[:,4],base_data[:,5]]
    Effective_Nuclear_Charge_Zeff = []
    Energy = base_data[:,1]
    Ionzation_Potential_D = []
    Ionzation_Potential_He3 = []
    Ionzation_Potential_He4 = []
    Ionzation_Potential_Impurities = [[] for _ in range(len(inp.IMPURITIES.keys()))]
    Ionzation_Potential_T = []
    Istopic_Mass_Number = []
    LRad_D = []
    LRad_He3 = []
    LRad_He4 = []
    LRad_Impurities = [[] for _ in range(len(inp.IMPURITIES.keys()))]
    LRad_T = []
    Power_Ionization_Potential_Total = []
    Power_Nuclear = []
    Power_AuxiliaryHeating = []
    Power_Boundary = []
    Power_Bremmstrahlung = []
    Power_Joule = []
    Power_Radiation = []
    Power_Reactions = [[] for v in inp.REACTIONS.values() if v]
    Power_Sink = []
    Power_Source = []
    Power_Surface = []
    Power_ThresholdLH = []
    Rate_FusionSigmaV = [[] for v in inp.REACTIONS.values() if v]
    Rate_Reactions = [[] for v in inp.REACTIONS.values() if v]
    Resistivity = []
    Tau_E = []
    Tau_HMode = []
    Tau_LMode = []
    Time = base_data[:,0]
    Temperature = []
    
    # Derivated quantities
    for dd in base_data:
        t = dd[0]
        values = dd[1:]

        s = vector(values)
                           
        Effective_Nuclear_Charge_Zeff.append(s.getZeff())
        Ionzation_Potential_D.append(s.ElementIonizationLosses('D'))
        Ionzation_Potential_He4.append(s.ElementIonizationLosses('He4'))
        Ionzation_Potential_T.append(s.ElementIonizationLosses('T'))
        Istopic_Mass_Number.append(s.isotopicMassNumber())
        LRad_D.append(s.LineRadiationPEC(ut.PEC['D']))
        LRad_He4.append(s.LineRadiationPEC(ut.PEC['He4']))
        LRad_T.append(s.LineRadiationPEC(ut.PEC['T']))
        Power_Ionization_Potential_Total.append(s.IonizationPotential())
        Power_Nuclear.append(sum(s.NuclearPower(reaction) if inp.REACTIONS[reaction] else 0. for reaction in inp.REACTIONS.keys()))
        Power_AuxiliaryHeating.append(s.AuxiliaryHeating(t))
        Power_Boundary.append(s.BoundaryPowerLosses(t))
        Power_Bremmstrahlung.append(s.Bremsstrahlung())
        Power_Joule.append(s.JoulePower())
        Power_Radiation.append(s.RadiatedPower())
        Power_Sink.append(s.PowerSink(t))
        Power_Source.append(s.PowerSource(t))
        Power_Surface.append(s.surfacePowerLoss(t))
        Power_ThresholdLH.append(s.getHLPowerThreshold())
        Resistivity.append(s.getResistivity())
        Tau_E.append(s.getTauE(t))
        Tau_HMode.append(s.HmodeConfinementTime(t))
        Tau_LMode.append(s.LmodeConfinementTime(t))
        Temperature.append(s.TeV)
        
        # Plotting He3 specific data
        if Density_He3:
            Ionzation_Potential_He3.append(s.ElementIonizationLosses('He3'))
            LRad_He3.append(s.LineRadiationPEC(ut.PEC['He3']))
        
        # Getting only the selected reactions data        
        index = 0
        for reaction in inp.REACTIONS.keys():
            if inp.REACTIONS[reaction]:
                Rate_FusionSigmaV[index].append(s.getReactionRate(reaction))
                Rate_Reactions[index].append(s.getReactionsSecond(reaction))
                Power_Reactions[index].append(s.NuclearPower(reaction))
                index += 1

        # Impurities
        for i in range(len(inp.IMPURITIES.keys())):
            element = list(inp.IMPURITIES.keys())[i]
            Ionzation_Potential_Impurities[i].append(s.ElementIonizationLosses(element))
            LRad_Impurities[i].append(s.LineRadiationPEC(ut.PEC[element]))
            
            # Adding also the total impurity density
            tot = 0
            for density in s.nImpurities[element]:
                tot += density
                
            Density_Impurities[i].append([*s.nImpurities[element], tot])
    
    # Labels for nuclear reactions
    labels_react = [k for k,v in inp.REACTIONS.items() if v]
    
    # Plasma total composition
    Plasma_Composition = [[sum(ion) for ion in zip(*element)] for element in [Density_D, Density_T, Density_He4]]
    Plasma_Composition_labels = ['D', 'T', 'He$_4$']
    
    # Line radiation
    LRad = [LRad_D, LRad_T, LRad_He4]
    
    # Ionization potential
    Ionization_Potential = [Ionzation_Potential_D, Ionzation_Potential_T, Ionzation_Potential_He4]
    
    # He3
    if Density_He3:
        Plasma_Composition.append([sum(ion) for ion in zip(*Density_He3)])
        Ionization_Potential.append(Ionzation_Potential_He3)
        LRad.append(LRad_He3)
        Plasma_Composition_labels.append('He$_3$')
        
    # Impurities
    for i in range(len(inp.IMPURITIES.keys())):
        element = list(inp.IMPURITIES.keys())[i]
        Plasma_Composition.append(np.transpose(Density_Impurities[i])[-1])
        Plasma_Composition_labels.append(element)
    
    # Definition of the plots
    plots = {
        'Temperature evolution': [Time,Temperature,'Temperature (eV)','',False],
        'Total energy content': [Time,Energy,'Plasma energy content (J)','',False],
        'Plasma composition': [Time,Plasma_Composition,'Particle density (m$^{-3}$)',Plasma_Composition_labels,True],
        'Energy balance details': [Time, [Power_Sink, Power_Source], 'Power (W)', ['Sinks','Sources'],False],
        'Energy sources detais': [Time, [Power_Nuclear, Power_AuxiliaryHeating, Power_Joule], 'Power (W)', ['Alpha power','Auxiliary heating','Joule power'],True],
        'Energy sinks details': [Time, [Power_Boundary, Power_Bremmstrahlung, Power_Radiation, Power_Ionization_Potential_Total], 'Power (W)', ['Boundary losses', 'Bremmstrahlung radiation', 'Line radiation', 'Ionization potential'],True],
        'Line radiation details': [Time, LRad+[LRad_Imp for LRad_Imp in LRad_Impurities], 'Power (W)', Plasma_Composition_labels,True],
        'Ionization potential losses': [Time, Ionization_Potential+[Ion_imp for Ion_imp in Ionzation_Potential_Impurities], 'Power (W)', Plasma_Composition_labels,True],
        'Plasma confinement mode': [Time, [Power_Surface, Power_ThresholdLH], 'Power(W)', ['Surface power loss','Threshold for H-mode'],False],
        'Plasma confinement time': [Time, [Tau_LMode, Tau_HMode, Tau_E], 'Confinement time (s)', [r'$\tau_{L(ITER89-P)}$',r'$\tau_{H(98y2)}$',r'$\tau_E$'],False],
        'Deuterium density evolution': [Time, Density_D, 'Particle density (m$^{-3}$)', ['D$^0$', 'D$^+$'],True],
        'Tritium density evolution': [Time, Density_T, 'Particle density (m$^{-3}$)', ['T$^0$', 'T$^+$'],True],
        '4-Helium density evolution': [Time, Density_He4, 'Particle density (m$^{-3}$)', ['He$_4^0$','He$_4^+$','He$_4^{2+}$'],True],
        r'Fusion rate $\left\langle \sigma v \right\rangle$': [Time, Rate_FusionSigmaV, r'$\left\langle \sigma v \right\rangle$ (m$^3$ s$^{-1}$)', labels_react,True],
        'Fusion reactions per second': [Time, Rate_Reactions, 'Reactions per second (s$^{-1}$)', labels_react,True],
        'Effective nuclear charge evolution': [Time, Effective_Nuclear_Charge_Zeff, r'Z$_{eff}$ (-)','',False],
        'Plasma resistivity evolution': [Time, Resistivity, r'$\rho$ ($\Omega$ m)','',True],
        'Average isotopic mass number': [Time, Istopic_Mass_Number, 'M [-]','',False]
    }            
    
    # Impurity details
    for i in range(len(inp.IMPURITIES.keys())):
        element = list(inp.IMPURITIES.keys())[i]
        labels = []
        for charge in range(ut.IMPURITIES_CHARGE_STATES[element]):
            labels.append(f'{element} {charge}' if charge > ut.IMPURITIES_CHARGE_STATES[element] - 7 else None)
        labels.append(f'{element} tot')
        plots[f'{element} density evolution'] = [Time, np.transpose(Density_Impurities[i]), 'Particle density (m$^{-3}$)',labels,False]
    
    # He3
    if Density_He3:
        plots['3-Helium density evolution'] = [Time, Density_He3, 'Particle density (m$^{-3}$)', ['He$_3^0$','He$_3^+$','He$_3^{2+}$'],True]    
        
    # Parassitic reactions
    if len(Power_Reactions) > 1:
        plots['Nuclear reactions details'] = [Time, Power_Reactions, 'Power (W)', labels_react, True]
    
    # Plots
    for i in range(len(plots.keys())):
        # Plotting every 10 plots to not overload the memory
        if i % 10 == 0:
            plt.plot()
        
        title = list(plots.keys())[i]    
        x_value, y_values, y_legend, y_labels, y_log = list(plots.values())[i]
        
        plt.figure(i)
        plt.title(title)
        plt.plot(x_value, np.transpose(y_values), label=y_labels)
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel(y_legend)
        plt.grid(True, 'both')
        if y_log:
            plt.semilogy()
            
    plt.plot()
