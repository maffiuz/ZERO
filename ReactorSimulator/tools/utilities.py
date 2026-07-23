# This file contains useful constant values, and coefficients datasets 
import os
import userDefined.user_inputs as inp
from tokamak.adas import reads_acd as reads_acd
from tokamak.adas import reads_acd as reads_scd
from tokamak.adas import reads_pec as reads_pec 
import tokamak.fusionrates as fr

# Useful constants
NSMALL = 1.0e-15 

CB = 5.35e3             # Bremmstrahlung radiation constant
KAPPA = 1.75               # -, plasma elongation parameter
TAU_SMOOTHING = 0.0001  # Time constant for smoothing constant injection rates of impurities

# Elements parameters, filled in run_initializer.py
IMPURITIES_CHARGE_STATES = {}
IMPURITIES_MASS_NUMBER = {}

# Parameter for forcing confinment time
TAU_FIXED = .2

# Parameters for handling the L to H mode switch
INTER_PREC = 5.e2       # -, velocity of interpolation for L-H mode switch
INTER_THRESH = 1.e-3    # -, relative error on the interpolating arctan function for L-H mode switch

# Referencing databases
# H
recHfile = os.path.join('adas','H','acd96_h.dat')
ionHfile = os.path.join('adas','H','scd96_h.dat')
# T
recTfile = os.path.join('adas','H','acd96_t.dat')
ionTfile = os.path.join('adas','H','scd96_t.dat')
# He
recHefile = os.path.join('adas','He','acd74_he.dat')
ionHefile = os.path.join('adas','He','scd74_he.dat')
# Ar
recArfile = os.path.join('adas','Ar','acd89_ar.dat')
ionArfile = os.path.join('adas','Ar','scd89_ar.dat')

# PEC
pec12_H_pfufile = os.path.join('adas','H','pec12#h_pju#h0.dat')
pec96_He_file = os.path.join('adas','He','pec96#he_pju#he0.dat')
pec96_He1_file = os.path.join('adas','He','pec96#he_pju#he1.dat')

def getPEC_data(element: str) -> list:
    # Setting appropriate path in the right folder (download has to be manual)
    path = f'adas/{element}'
    
    # Fetching all 'pec... .dat'
    pec_data_files = [
        f'{path}/{f}' for f in os.listdir(path)
        if os.path.isfile(os.path.join(path, f))
        and 'pec' in f
        and f.endswith('.dat')
    ]   

    pec_data = list(range(len(pec_data_files)))

    # Sorting pec files based on charge state
    for f in pec_data_files:
        # Handling two digits
        try:
            charge = int(f[-6]+f[-5])
        except:
            charge = int(f[-5])

        # Creating readers
        pec_data[charge] = reads_pec(f)

    return pec_data


# Building readers
# H
recH = reads_acd(recHfile)
ionH = reads_scd(ionHfile)
# T
recT = reads_acd(recTfile)
ionT = reads_scd(ionTfile)
# He
recHe = reads_acd(recHefile)
ionHe = reads_scd(ionHefile)
# Ar
recAr = reads_acd(recArfile)
ionAr = reads_scd(ionArfile)

# PEC
pec12_H_pju = reads_pec(pec12_H_pfufile)
pec96_He = reads_pec(pec96_He_file)
pec96_He1 = reads_pec(pec96_He1_file)


# Dictionaries constructed to properly access the recombination/ionization functions
REC = {
    'D':[recH,lambda s,Z:s.nD[Z]],
    'T':[recT,lambda s,Z:s.nT[Z]],
    'He4':[recHe,lambda s,Z:s.nHe4[Z]],
    'He3':[recHe,lambda s,Z:s.nHe3[Z]],
    'Ar':[recAr,lambda s,Z:s.nImpurities['Ar'][Z]]
}

ION = {
    'D':[ionH,lambda s,Z:s.nD[Z-1]],
    'T':[ionT,lambda s,Z:s.nT[Z-1]],
    'He4':[ionHe,lambda s,Z:s.nHe4[Z-1]],
    'He3':[ionHe,lambda s,Z:s.nHe3[Z-1]],
    'Ar':[ionAr,lambda s,Z:s.nImpurities['Ar'][Z-1]]
}

# Dictionary constructed to properly access the refueling function
FUEL_IN = {
    'D':inp.D0_in,
    'T':inp.T0_in
}

# Dictionary constructed to properly access the ion mass loss function
MASS_LOSS_TAU = {
    'D':lambda s,Z:s.nD[Z]/5.0,
    'T':lambda s,Z:s.nT[Z]/5.0,
    'He4':lambda s,Z:s.nHe4[Z]/5.0,
    'He3':lambda s,Z:s.nHe3[Z]/5.0,
    'Imp':lambda el,s,Z:s.nImpurities[el][Z]/5.0
}

# Dictionary constructed to properly access the recimbined particles return function
# % of lost particles
BOUNCE = {
    'D':0.99,
    'T':0.99,
    'He4':0.15,
    'He3':0.15,
    'Imp':0.99
}

# Dictionary constructed to properly access the line radiation (PEC) function
PEC = {
    'D':[getPEC_data('H'), lambda s,Z:s.nD[Z]],
    'T':[getPEC_data('H'), lambda s,Z:s.nT[Z]],
    'He4':[getPEC_data('He'), lambda s,Z:s.nHe4[Z]],
    'He3':[getPEC_data('He'), lambda s,Z:s.nHe3[Z]],
    'Ar':[getPEC_data('Ar'), lambda s,Z:s.nImpurities['Ar'][Z]],
}

# Dictionary to properly access reaction rates
REACTION_RATES = {
    'H3(D,N)HE4': lambda TeV: fr.rate_H3_DN_HE4(TeV),
    'H2(D,N)HE3': lambda TeV: fr.rate_H2_DN_HE3(TeV),
    'H2(D,P)H3': lambda TeV: fr.rate_H2_DP_H3(TeV),
    'HE3(D,P)HE4': lambda TeV: fr.rate_HE3_DP_HE4(TeV)
}

# Dictionary to properly access ions for fusion reactions
REACTION_IONS = {
    'H3(D,N)HE4': lambda s: [s.nD[1], s.nT[1]],
    'H2(D,N)HE3': lambda s: [s.nD[1], s.nD[1]],
    'H2(D,P)H3': lambda s: [s.nD[1], s.nD[1]],
    'HE3(D,P)HE4': lambda s: [s.nHe3[2], s.nD[1]]
}

# Dictionary to properly access power yields for fusion reactions [MeV]
REACTION_POWER = {
    'H3(D,N)HE4': 17.589 * 0.2,
    'H2(D,N)HE3': 3.269 * 0.25,
    'H2(D,P)H3': 4.033,
    'HE3(D,P)HE4': 18.353
}

# Ionization potentials
with open('data/ionization_potentials','r') as f:
    lines = f.readlines()

POTENTIALS = {
    'D':lines[0].split(),
    'T':lines[0].split(),
    'He4':lines[1].split(),
    'He3':lines[1].split(),
}
