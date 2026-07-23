# %%
import scipy.constants as cs
from numpy import exp as exp

def TeV2T9(TeV):
    # Convert temperature in eV to local units (1e9K)
    T9 = 1e-9*TeV*cs.eV/cs.k
    return T9

def T913(T9):
    return pow(T9,1./3.)

def T923(T9):
    return pow(T9,2./3.)

def T932(T9):
    return pow(T9,3./2.)

def T943(T9):
    return pow(T9,4./3.)

def T953(T9):
    return pow(T9,5./3.)

def rate_H3_DN_HE4(TeV: float) -> float:
    
    # Return the rate coefficient for the H3(D,N)HE4 reaction
    # D + T -> He4 + n
    T9 = TeV2T9(TeV)
    T9_23 = T923(T9)
    T9_13 = T913(T9)
    T9_43 = T943(T9)
    T9_53 = T953(T9)
    # Convert the temperature to local units (1e9K) and evaluate relevant T powers.
     
    rr = 8.09e10/T9_23 * exp(-4.524/T9_13 - pow(T9/0.120,2.0)) * \
        (1.0 + 0.092*T9_13 + 1.80*T9_23 + 1.16*T9 + 10.52*T9_43 + 17.24*T9_53) + \
        8.73e8/T9_23*exp(-0.523/T9)
    # Correlation source highlighted in manual    
    
    rr_cgs = rr/cs.N_A
    rr_mks = rr_cgs * 1.e-6
    # Convert to mks units
    
    return rr_mks

def rate_H2_DP_H3(TeV: float) -> float:
        
    # Return the rate coefficient for the H2(D,P)H3 reaction
    # D + D -> T + p
    T9 = TeV2T9(TeV)
    T9_23 = T923(T9)
    T9_13 = T913(T9)
    T9_43 = T943(T9)
    T9_53 = T953(T9)
    # Convert the temperature to local units (1e9K) and evaluate relevant T powers.
     
    rr = 4.13e8/T9_23 * exp(-4.258/T9_13) * \
        (1.0 + 0.098*T9_13 + 4.39e-2*T9_23 + 3.01e-2*T9 + 0.543*T9_43 + 0.946*T9_53) + \
        1.73*exp(-46.798/T9)
    # Correlation source highlighted in manual    
    
    rr_cgs = rr/cs.N_A
    rr_mks = rr_cgs * 1.e-6
    # Convert to mks units
    
    return rr_mks

def rate_H2_DN_HE3(TeV: float) -> float:
        
    # Return the rate coefficient for the H2(D,N)HE3 reaction
    # D + D -> He3 + n
    T9 = TeV2T9(TeV)
    T9_23 = T923(T9)
    T9_13 = T913(T9)
    T9_43 = T943(T9)
    T9_53 = T953(T9)
    # Convert the temperature to local units (1e9K) and evaluate relevant T powers.
     
    rr = 3.88e8/T9_23 * exp(-4.258/T9_13) * \
        (1.0 + 0.098*T9_13 + 0.418*T9_23 + 0.287*T9 + 0.638*T9_43 + 1.112*T9_53) + \
        1.73*exp(-37.935/T9)
    # Correlation source highlighted in manual    
    
    rr_cgs = rr/cs.N_A
    rr_mks = rr_cgs * 1.e-6
    # Convert to mks units
    
    return rr_mks

def rate_HE3_DP_HE4(TeV: float) -> float:
        
    # Return the rate coefficient for the HE3(D,P)HE4 reaction
    # He3 + D -> He4 + p
    T9 = TeV2T9(TeV)
    T9_23 = T923(T9)
    T9_13 = T913(T9)
    T9_43 = T943(T9)
    T9_53 = T953(T9)
    # Convert the temperature to local units (1e9K) and evaluate relevant T powers.
     
    rr = 5.86e10/T9_23 * exp(-7.181/T9_13 - (T9/0.315)**2) * \
        (1.0 + 0.058*T9_13 + 0.142*T9_23 + 5.78e-2*T9 + 2.25*T9_43 + 2.32*T9_53) + \
        4.36e8/pow(T9, 1/2) * exp(-1.720/T9)
    # Correlation source highlighted in manual    
    
    rr_cgs = rr/cs.N_A
    rr_mks = rr_cgs * 1.e-6
    # Convert to mks units
    
    return rr_mks
    

def main():
    # Useful to visualize the rate dependency from T
    import matplotlib.pyplot as plt
    import numpy as np
    
    TT = np.logspace(2,5,500)
    
    rates = [[rate_H3_DN_HE4(TeV), rate_H2_DP_H3(TeV), rate_H2_DN_HE3(TeV), rate_HE3_DP_HE4(TeV)] for TeV in TT]
    labels = ['D+T->He4','D+D->T','D+D->He3','He3+D->He4']
    Tk = [TeV2T9(TeV)*1e9 for TeV in TT ]
    
    fg = plt.figure('Reaction rates')
    ax = plt.gca()
    plt.plot(Tk,rates,label=labels)
    plt.grid()
    plt.legend()
    ax.set_xlabel('T (K)')
    ax.set_ylabel('rate (m3 s-1)')
    ax.set_xscale('log')
    ax.set_yscale('symlog', linthresh=1e-25)
    ax.set_yticks([0, 1e-25, 1e-24, 1e-23, 1e-22, 1e-21, 1e-20],['0', '1e-25', '1e-24', '1e-23', '1e-22', '1e-21', '1e-20'])
    ax.set_ylim(bottom=0)
    plt.show()    
    
if __name__ == "__main__":
    main()
        
# %%
