# This file contains all the methods to create the output reports of the script
import tools.utilities as ut
import userDefined.user_inputs as inp
import datetime

# Writer function     
def writer(t,y):
    # Main output file
    with open('out/simulator_data','w') as f:
        # Check if it is needed to add He3
        if inp.REACTIONS['H2(D,N)HE3'] or inp.REACTIONS['HE3(D,P)HE4']:
            index = 11
            f.write('   time (s)      Etot (j)      nD (m-3)      nD+ (m-3)     nT (m-3)     nT+ (m-3)     nHe4 (m-3)    nHe4+ (m-3)  nHe4++ (m-3)   nHe3 (m-3)    nHe3+ (m-3)  nHe3++ (m-3)\n'+'='*14*12+'\n')
            
            for k in range(0, len(t)):
                f.write(f'{t[k]:13.6e} ')       # Time         (s)
                f.write(f'{y[0][k]:13.6e} ')    # Total energy (J)
                f.write(f'{y[1][k]:13.6e} ')    # Neutral D density (m-3)
                f.write(f'{y[2][k]:13.6e} ')    # Ionized D density (m-3)
                f.write(f'{y[3][k]:13.6e} ')    # Neutral T density (m-3)
                f.write(f'{y[4][k]:13.6e} ')    # Ionized T density (m-3)
                f.write(f'{y[5][k]:13.6e} ')    # Neutral He4 density (m-3)
                f.write(f'{y[6][k]:13.6e} ')    # He4 + density (m-3)
                f.write(f'{y[7][k]:13.6e} ')    # He4 2+ density (m-3)
                f.write(f'{y[8][k]:13.6e} ')    # Neutral He3 density (m-3)
                f.write(f'{y[9][k]:13.6e} ')    # He3 + density (m-3)
                f.write(f'{y[10][k]:13.6e} ')   # He3 2+ density (m-3)
                f.write('\n')
                
        else:
            index = 8
            f.write('   time (s)      Etot (j)      nD (m-3)      nD+ (m-3)     nT (m-3)     nT+ (m-3)     nHe4 (m-3)    nHe4+ (m-3)  nHe4++ (m-3)\n'+'='*14*9+'\n')
            
            for k in range(0, len(t)):
                f.write(f'{t[k]:13.6e} ')       # Time         (s)
                f.write(f'{y[0][k]:13.6e} ')    # Total energy (J)
                f.write(f'{y[1][k]:13.6e} ')    # Neutral D density (m-3)
                f.write(f'{y[2][k]:13.6e} ')    # Ionized D density (m-3)
                f.write(f'{y[3][k]:13.6e} ')    # Neutral T density (m-3)
                f.write(f'{y[4][k]:13.6e} ')    # Ionized T density (m-3)
                f.write(f'{y[5][k]:13.6e} ')    # Neutral He4 density (m-3)
                f.write(f'{y[6][k]:13.6e} ')    # He4 + density (m-3)
                f.write(f'{y[7][k]:13.6e} ')    # He4 2+ density (m-3)
                f.write('\n')
    
    for element in inp.IMPURITIES.keys():
        index = impuritiesWriter(t,y,element,index)
            
# Specific function to output impurities data
def impuritiesWriter(t,y, element: str, index: int) -> int:
    states = ut.IMPURITIES_CHARGE_STATES[element]
    path = f'out/impurities_{element}_data'
    # Impurities
    with open(path,'w') as f:
        f.write(f'{element} impurity densities (m-3), all charge states ({states})\n')
        f.write('='*14*(states+1)+'\n')
        f.write('   time (s)   ' + ''.join([f'    Z = {z:^2}    ' for z in range(states)]) + '\n')
        f.write('='*14*(states+1)+'\n')
               
        for k in range(0, len(t)):
            f.write(f'{t[k]:13.6e} ')     # Time         (s)
            
            for j in range(index, index+states):
                f.write(f'{y[j][k]:13.6e} ')  # Density (m-3)
            
            f.write('\n')
            
    return index + states

# Writing a report with all input data. This can be itself used as input
def startupReport():
    current_time = datetime.datetime.now() 
    
    with open(f'out/run_report_{inp.RUN_NAME}','w', encoding="utf-8") as f:           
        f.write(f'[{str(current_time)}]\n')
        f.write('='*72+'\n'+'--- RUN REPORT '+57*'-'+'\n'+'='*72+'\n')
        f.write(
            f"\n--- RUN NAME:\n\
        {inp.RUN_NAME}\n\n\n\
--- TOKAMAK CHAMBER\n\n\
        R = {inp.R:<16.2f}[m]    major plasma radius\n\
        a = {inp.a:<16.2f}[m]    minor plasma radius\n\n\n\
--- TIME FRAME\n\n\
        t_start = {inp.t_start:<10.2f}[s]    simulation starting time (this is just for reference)\n\
        t_end = {inp.t_end:<12.2f}[s]    simulation ending time\n\
        n_points = {inp.n_points:<9}[-]    number of points on which the solver projects the solution\n\n\n\
--- OPERATIONAL PARAMETERS\n\n\
    ------- Auxiliary power:\n\
        W_ext = {inp.Aux_heat_operat:<12.2e}[W]    auxiliary heating\n\n\
    ------- Startup:\n\
        W_ext = {inp.Aux_heat_start:<12.2e}[W]    auxiliary heating during startup\n\
        t_down = {inp.TimeDown:<11.2f}[s]    startup time\n\n\
    ------- Confinment:\n\
        BT = {inp.BT:<15.2f}[T]    toroidal magnetic field component\n\
        IP = {inp.IP:<15.2e}[A]    induced current in the plasma\n\n\
    ------- Refueling:\n\
        t_ref = {inp.t_ref:<12.2f}[s]    refueling starting time\n\
        D0_in = {inp.D0_in:<12.2e}[s-1]  deuterium refueling rate\n\
        T0_in = {inp.T0_in:<12.2e}[s-1]  tritium refueling rate\n\n\n\
--- STARTING CONDITIONS\n\n\
        E = {inp.E_START:<16.2e}[J]    total initial energy content\n\
        nD0 = {inp.ND0_START:<14.2e}[m-3]  initial neutral deuterium particle density\n\
        nD1 = {inp.ND1_START:<14.2e}[m-3]  initial ionized deuterium particle density\n\
        nT0 = {inp.NT0_START:<14.2e}[m-3]  initial neutral tritium particle density\n\
        nT1 = {inp.NT1_START:<14.2e}[m-3]  initial ionized tritium particle density\n\
        nHe4_0 = {inp.NHE4_START:<11.2e}[m-3]  initial neutral He4 particle density\n\
        nHe4_0 = {inp.NHE41_START:<11.2e}[m-3]  initial ionized He4 particle density\n\
        nHe4_0 = {inp.NHE42_START:<11.2e}[m-3]  initial alpha particle density\n\n\
"
        )
        
        # Reactions 
        f.write(f'\n--- REACTIONS\n\n\
        H3(D,N)HE4 = {inp.REACTIONS['H3(D,N)HE4']}\n\
        H2(D,N)HE3 = {inp.REACTIONS['H2(D,N)HE3']}\n\
        H2(D,P)H3 = {inp.REACTIONS['H2(D,P)H3']}\n\
        HE3(D,P)HE4 = {inp.REACTIONS['HE3(D,P)HE4']}\n\n')
        
        # Handling impurities
        if inp.IMPURITIES:
            f.write('\n--- IMPURITIES\n\n')
        
        for el in inp.IMPURITIES.keys():
            f.write(f"\
    Impurity element: {el}\n\
        n{el} = {inp.IMPURITIES[el][0]:<14.2e}[m-3]  initial neutral {el} particle density\n\
        inj_rate = {inp.IMPURITIES[el][1]:<9.2e}[s-1]  {el} injection rate\n\
        inj_start = {inp.IMPURITIES[el][2]:<8.2f}[s]    {el} injection starting time\n\
        inj_dur = {inp.IMPURITIES[el][3]:<10.2f}[s]    {el} injection duration\n\n")
            
        # Bottom line
        f.write('\n'+'='*72)
        f.write('\nTHIS LINE MARKS THE END OF THE REPORT. ANY COMMENTS CAN BE WRITTEN BELOW\n')
        f.write("\nNote: don't change the format of the report file if you wish to use it as \ninput for future runs\n")
