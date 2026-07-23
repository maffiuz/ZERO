# Collection of variables that users can modify
import scipy.constants as cs

# Run name
RUN_NAME = 'Test'

# Time frame
t_start = 0. 
t_end = 300
n_points = 1000

# Reactor chamber dimentions
a = 2.                 # m, plasma section radius
R = 6.2                # m, tokamak major radius

# Operational variables
t_ref = 60              # s, refueling starting time
D0_in = 5.e20           # s-1, deuterium particle refueling rate 
T0_in = 5.e20           # s-1, tritium particle refueling rate
TimeDown = 60          # s, time of startup auxiliary heating shut down
Aux_heat_start = 200.e6   # W, auxiliary heating during startup
Aux_heat_operat = 10.e6   # W, auxiliary heating during operation
IP = 12.5e6              # A, induced current in the plasma
BT = 11.8               # T, toroidal field component

# Starting conditions
E_START = 1.e-2*cs.eV*1.e20 # J, total initial energy content
ND0_START = 1.0e20      # m-3, initial neutral deuterium particle density
ND1_START = 1.          # m-3, initial ionized deuterium particle density
NT0_START = 1.0e20      # m-3, initial neutral tritium particle density
NT1_START = 1.          # m-3, initial ionized tritium particle density
NHE4_START = 0.         # m-3, initial neutral He4 particle density
NHE41_START = 0.        # m-3, initial ionized He4 particle density
NHE42_START = 0.        # m-3, initial alpha particle density

# Impurities
IMPURITIES = {          # starting density, injection rate, injection start time, injection time 
    'Ar': [0, 4.7e17, 100, 10],
}

# Reactions
REACTIONS = {
    'H3(D,N)HE4': True,
    'H2(D,N)HE3': True,   
    'H2(D,P)H3': True,     
    'HE3(D,P)HE4': True        
}
