# This file contains the vector class, used to access and manipulate data
# This class is called every step of the iteration as it allows the computation
# of key tokamak parameters, given the parameteres contained in the main equations of the solver
import scipy.constants as cs
import numpy as np
import userDefined.user_inputs as inp
import tools.utilities as ut

# Class used to read and operate on the unknowns vector
class vector:
    # Constructor
    def __init__(self, s: list[float]):
        self.Etot = s[0]                                                    # Total energy content
        
        # Ion densities are stored in tuples (no need to change them once the class is instanced)
        # To access the charghe state simply use the index (i.e. charge state Z=0 -> index 0)
        self.nD = (s[1], s[2])                                              # D ion densities
        self.nT = (s[3], s[4])                                              # T ion densities
        self.nHe4 = tuple(s[5:8])                                           # He4 ion densities
        
        # End of base variables
        index = 8
        self.nHe3 = ()
        # Activate He3 only in presence of parassitic reactions that may involve it
        if inp.REACTIONS['H2(D,N)HE3'] or inp.REACTIONS['HE3(D,P)HE4']:
            self.nHe3 = tuple(s[index:index+3])                             # He3 ion densities
            index += 3
            
        # Impurities
        self.nImpurities = {}
        for element in inp.IMPURITIES.keys():
            self.nImpurities[element] = tuple(s[index:index+ut.IMPURITIES_CHARGE_STATES[element]])
            index += ut.IMPURITIES_CHARGE_STATES[element]
                    
        # Auxiliary useful parameters
        self.ne = self.setne()                                              # Electron density
        self.V = (2*cs.pi*inp.R) * (pow(inp.a,2.0)*cs.pi)                   # Plasma volume
        self.TeV = self.setTeV()                                            # Plasma temperature (eV)
    
    # Function to set the temperature in eV
    def setTeV(self) -> float:
        TJ = 2/3 * self.Etot / self.V / (self.getTotaln())
        TeV = TJ/cs.eV
        # Just to make sure that the T never actually goes to 0
        # Avoids overflow error with log functions
        return max(TeV,ut.NSMALL)  
    
    # Function to set the number of electrons
    def setne(self) -> float:
        ne = self.nD[1] + self.nT[1] + self.nHe4[1] + 2*self.nHe4[2]   # Base electron density
        
        if self.nHe3:
            ne += self.nHe3[1] + 2*self.nHe3[2]
            
        # Electrons from impurity ionization
        for impurity in self.nImpurities.values():
            for i in range(len(impurity)):
                ne += i*impurity[i]
            
        return ne
        
    # Function to get the log10 of the temperature in eV    
    def getTeVlog10(self) -> float:
        return np.log10(self.TeV)
    
    # Total particle density [m-3]
    def getTotaln(self) -> float:
        base_sum = self.nD[0] + self.nD[1] + self.nT[0] + self.nD[1] + self.nHe4[0] + self.nHe4[1] + self.nHe4[2] + self.ne
        if self.nHe3:
            base_sum += sum(self.nHe3)
        impurity_sum = sum(sum(impurity) for impurity in self.nImpurities.values())
        return base_sum + impurity_sum
    
    # Function to get the reaction rate for the given Temperature [m3 s-1]   
    def getReactionRate(self, reaction: str) -> float:
        return ut.REACTION_RATES[reaction](self.TeV) if inp.REACTIONS[reaction] else 0.
    
    # Functions to get the reactions per second, given the rate parameter [m-3 s-1]   
    def getReactionsSecond(self, reaction: str) -> float:
        element1, element2 = ut.REACTION_IONS[reaction](self)
        rate = self.getReactionRate(reaction)
        return element1 * element2 * rate 
    
    # Generic recombination/ionization function [m-3 s-1]
    def RecIonHandler(self,x: list, Z: int) -> float:
        """Generic recombination/ionization function

        Args:
            x (list): the appropriate reader for the function: REC/ION with specified key the element
            Z (int): charge state. Note that Z is the START for recombinations and is the TARGET for ionization

        Returns:
            float: the recombination/ionization mass rate [m-3 s-1]
        """
        # Get the appropriate data file, charge state and density
        data = x[0]
        density = x[1](self,Z)
        charge = f'Z={Z}'
        
        Tlog10 = self.getTeVlog10()
        nlog10 = np.log10(max(density,ut.NSMALL))
        rate = pow(10,data['EVALUATE'][charge].ev(Tlog10,nlog10))
        
        # Return the recombination/ionization rate
        ion = rate * self.ne * density
        return ion
     
    # Generic refueling fucntion
    def Refueling(self,t: float,el: str) -> float:
        '''t_start = inp.t_ref
        rate = ut.FUEL_IN[el]/self.V
        
        # A sudden jump in change rate causes numerical instability, a smooth transition is necessary
        if t<t_start:
            return 0.
        elif t<t_start + ut.TAU_SMOOTHING:
            smoothness = (t-t_start)/ut.TAU_SMOOTHING
        else: 
            smoothness = 1
             
        return smoothness*rate'''
        
        return ut.FUEL_IN[el]/self.V if t>inp.t_ref else 0.
    
    # Generic mass loss function
    def MassLoss(self,el: str, Z: int, specie = None) -> float:
        if el == 'Imp':
            return ut.MASS_LOSS_TAU['Imp'](specie,self,Z)
        else:
            return ut.MASS_LOSS_TAU[el](self,Z)
    
    # Generic recombined return function
    def ParticleRecycling(self,el: str, Z: int, specie = None) -> float:
        if el == 'Imp':
            return ut.BOUNCE[el]*self.MassLoss(el,Z,specie)
        else:
            return ut.BOUNCE[el]*self.MassLoss(el,Z)
       
    # Heating power in the plasma produced by nuclear fusion reactions [W]
    def NuclearPower(self, reaction: str) -> float:
        reactions_second = self.getReactionsSecond(reaction)
        QJ_reaction = ut.REACTION_POWER[reaction] * (1.e6*cs.eV) * reactions_second * self.V
        return QJ_reaction
    
    # A completely fake value, to mimic external heating for a well defined time period [W]
    def AuxiliaryHeating(self,t: float) -> float:        
        return inp.Aux_heat_start if t <= inp.TimeDown else inp.Aux_heat_operat
    
    # Positive heating power in the plasma
    def PowerSource(self,t: float) -> float:
        IntrinsicHeating = sum(self.NuclearPower(reaction) if value else 0. for reaction,value in inp.REACTIONS.items())
        ExternalHeating = self.AuxiliaryHeating(t) 
        JouleEffect = self.JoulePower()
        
        Source = ExternalHeating + IntrinsicHeating + JouleEffect # W
        return Source

    # Advective contribution to the power loss [W]
    def BoundaryPowerLosses(self, t: float) -> float:
        return self.Etot/self.getTauE(t)
    
    # Generic PEC function for single charge state
    def ChargeStatePEC(self, data: list, density_fun: function, charge: int) -> float:
        density = density_fun(self,charge)        
        nlog10 = np.log10(max(density,ut.NSMALL))

        rateWrad = pow(10., data[charge]['REC'](self.getTeVlog10(),nlog10).item())  # J m3 s-1
        rad = rateWrad*self.ne*density*self.V
        return rad
    
    # Generic line radiation function
    def LineRadiationPEC(self, x: list) -> float:
        rad = 0
        data = x[0]
        density_fun = x[1]
        
        for z in range(len(data)):  
            rad += self.ChargeStatePEC(data, density_fun, z)
        return rad
    
    # Power losses due to radiation 
    def RadiatedPower(self) -> float:
        basePEC = self.LineRadiationPEC(ut.PEC['D']) +\
            self.LineRadiationPEC(ut.PEC['T']) +\
            self.LineRadiationPEC(ut.PEC['He4'])
        
        if self.nHe3:
            basePEC += self.LineRadiationPEC(ut.PEC['He3'])
        
        impuritiesPEC = 0
        for element in inp.IMPURITIES.keys():
            impuritiesPEC += self.LineRadiationPEC(ut.PEC[element])
                
        return basePEC + impuritiesPEC
    
    # Generic ionization losses function for single charge state
    def ChargeStateIonizationPotential(self, potential: float, ion_rate: float) -> float:
        return potential*ion_rate*self.V  # W
        
    # Generic ionization losses function per element 
    def ElementIonizationLosses(self, element: str) -> float:
        element_potentials = ut.POTENTIALS[element]
        loss = 0
        
        for z in range(len(element_potentials)):
            potential = float(element_potentials[z])*cs.eV             # J             
            ion_rate = self.RecIonHandler(ut.ION[element],z+1)    
            loss += self.ChargeStateIonizationPotential(potential, ion_rate)
            
        return loss
    
    # Power losses due to element ionization
    def IonizationPotential(self):
        loss = self.ElementIonizationLosses('D') +\
            self.ElementIonizationLosses('T') +\
            self.ElementIonizationLosses('He4')
            
        if self.nHe3:
            loss += self.ElementIonizationLosses('He3')
            
        for element in inp.IMPURITIES.keys():
            loss += self.ElementIonizationLosses(element)
            
        return loss
        
    # Negative heating power in the plasma
    def PowerSink(self, t: float) -> float:
        RadiationLosses = self.RadiatedPower()
        AdvectionLosses = self.BoundaryPowerLosses(t)
        BremsstrahlungLosses = self.Bremsstrahlung()
        IonizationLosses = self.IonizationPotential()
        
        Sink = RadiationLosses + AdvectionLosses + BremsstrahlungLosses + IonizationLosses #W
        return Sink
    
    # Computing the effective nuclear charge
    def getZeff(self) -> float:
        charge = self.nD[1] + self.nT[1] + self.nHe4[1] + 4*self.nHe4[2]
        
        if self.nHe3:
            charge += self.nHe3[1] + 4*self.nHe3[2]
        
        for element in self.nImpurities.values():
            for i in range(len(element)):
                charge += i*i*element[i]
                
        return charge/self.ne
    
    # Bremsstrahlung heat losses (ref. Freidberg)  [W]
    def Bremsstrahlung(self) -> float:
        loss = ut.CB*self.getZeff()*(self.ne/1e20)**2*pow(self.TeV/1000,1/2)
        return loss*self.V
    
    # Plasma resistivity (10.45 Friedberg)
    def getResistivity(self) -> float:
        eta = 6.5e-8 / pow((self.TeV/1000),3/2)
        return eta
    
    # Plasma Resistance
    def getResistance(self) -> float:
        Res = self.getResistivity()*2*cs.pi*inp.R/(cs.pi*inp.a**2)
        return Res 
    
    # Joule effect power
    def JoulePower(self) -> float:
        return self.getResistance()*inp.IP**2
    
    # Power threshold L - H mode (eq 3 - Martin)
    def getHLPowerThreshold(self) -> float:
        P = 2.15 * pow(self.ne * (1.e-20), 0.782) * pow(inp.BT, 0.772) * \
            pow(inp.a, 0.975) * pow(inp.R, 0.999) # MW
        return P * 1.e6
    
    # Computing the energy confinement time
    def getTauE(self,t: float) -> float:   
        # Interpolating L-mode and H-mode using normalized arctan
        # The interpolation starts at a defined relative error from PTreshold      
        
        surfP = self.surfacePowerLoss(t)
        P_threshold = self.getHLPowerThreshold()
        
        err_rel_threshold = (surfP-P_threshold)/P_threshold   # Relative error between the two        
        
        # Interpolating arctan function, going from -1/2 to 1/2
        interp_fun = np.arctan(ut.INTER_PREC*err_rel_threshold)/np.pi
        
        # Evaluation done on the arctan function
        # L-mode
        if interp_fun/0.5 < -(1 - ut.INTER_THRESH):
            tauE = self.LmodeConfinementTime(t)
        
        # Interpolation
        elif interp_fun/0.5 < 1 - ut.INTER_THRESH:
            tauE = (1/2 - interp_fun)*self.LmodeConfinementTime(t) +\
                (interp_fun + 1/2)*self.HmodeConfinementTime(t)
                
        # H-mode
        else:
            tauE = self.HmodeConfinementTime(t)
            
        return tauE if tauE > ut.TAU_FIXED else ut.TAU_FIXED
            
    # Surface net power loss
    def surfacePowerLoss(self, t: float) -> float:        
        return self.PowerSource(t) - self.Bremsstrahlung() - self.RadiatedPower()
        
    # Confinement time for H-Mode plasma
    def HmodeConfinementTime(self, t: float) -> float:
        # IPB98(y,2)
        HmodeConst = 0.0562 * pow(inp.IP * 1e-6, 0.93) * pow(inp.BT, 0.15) * pow(inp.R, 1.97) * pow(inp.a/inp.R, 0.58) * pow(ut.KAPPA, 0.78)  
        tauE = HmodeConst * pow(self.isotopicMassNumber(), 0.19) * pow(self.ne * 1.e-19, 0.41) * pow(self.surfacePowerLoss(t) * 1.e-6, -0.69)
                
        return tauE 
    
    # Confinement time for L-Mode plasma
    def LmodeConfinementTime(self, t: float) -> float:
        #ITER89-P
        LmodeConst = 0.048 * pow(inp.IP * 1e-6, 0.85) * pow(inp.R, 1.2) * pow(inp.a, 0.3) * pow(ut.KAPPA, 0.5) * pow(inp.BT, 0.2)
        tauE = LmodeConst * pow(self.isotopicMassNumber(), 0.5) * pow(self.ne * 1.e-20, 0.1) * pow(self.surfacePowerLoss(t) * 1.e-6, -0.5) 
             
        return tauE
    
    # Function for impurity injection
    def impurityInjection(self, t: float, impurity: str) -> float:
        rate = inp.IMPURITIES[impurity][1]/self.V
        t_start = inp.IMPURITIES[impurity][2]
        t_end = t_start + inp.IMPURITIES[impurity][3]
        
        # A sudden jump in change rate causes numerical instability, a smooth transition is necessary
        if t<t_start or t>t_end:
            return 0.
        elif t<t_start + ut.TAU_SMOOTHING:
            smoothness = (t-t_start)/ut.TAU_SMOOTHING
        elif t>t_end-ut.TAU_SMOOTHING:
            smoothness = (t_end-t)/ut.TAU_SMOOTHING
        else: 
            smoothness = 1
             
        return smoothness*rate
    
    # Average ion mass number
    def isotopicMassNumber(self) -> float:
        totalmass = sum(self.nD)*2 + sum(self.nT)*3 + sum(self.nHe4)*4
        
        if self.nHe3:
            totalmass += sum(self.nHe3)*3
               
        for element in self.nImpurities.keys():
            totalmass += sum(self.nImpurities[element])*ut.IMPURITIES_MASS_NUMBER[element]
        
        return totalmass/(self.getTotaln() - self.ne)
    