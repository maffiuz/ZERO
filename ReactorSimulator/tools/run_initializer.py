# This file generates starting conditions for the problem and sets user data
# either from user_inputs.py or from a previous run report
import os
from periodictable import elements
import userDefined.user_inputs as inp
import tools.utilities as ut

def run_initializer() -> list:   
    use_report = False
    input_report = ''
    
    # Fetching a report file in the main directory
    input_report_list = [
        f for f in os.listdir('.')
        if os.path.isfile(f'./{f}')
        and 'run_report_' in f
    ] 

    # Only if a report is present in the folder
    if input_report_list:
        # Asking user whether to use a report as input file
        while True:
            try:
                input_choice = input('Wish to use a report file as input file? [y/n]: ').lower()
                if input_choice not in ['y','n']:
                    print('Invalid choice')
                    continue
                elif input_choice == 'y':
                    use_report = True    
                break
                
            except Exception as e:
                print(f'Exception occurred: {e}')

        # If there is more then one report
        if use_report and len(input_report_list) > 1:
            for i in range(len(input_report_list)):
                print(f'{i} - {input_report_list[i]}')
                
            try:    
                input_report = input_report_list[int(input('Which report do you wish to use? (Insert the number only): '))]  
            except Exception as e:
                print(f'Exception occurred: {e}')

        else:
            input_report = input_report_list[0]
    
    # Reading the report and overwriting user_inputs.py 
    s0 = []     
    
    if use_report:
        print(f'{input_report} will be used as input file...') 
        s0 = read_report(input_report)
    else:
        print('Using the values contained in userDefined/user_inputs.py ...')
        # Starting conditions          
        s0.append(inp.E_START)          # Total energy content (J)
        s0.append(inp.ND0_START)        # Deuterium neutral density
        s0.append(inp.ND1_START)        # Deuterium ion density
        s0.append(inp.NT0_START)        # Tritium neutral density
        s0.append(inp.NT1_START)        # Tritium ion density
        s0.append(inp.NHE4_START)       # He4 neutral density
        s0.append(inp.NHE41_START)      # He4 + ion density
        s0.append(inp.NHE42_START)      # He4 2+ ion density
        
    # Add He3 if needed
    if inp.REACTIONS['H2(D,N)HE3'] or inp.REACTIONS['HE3(D,P)HE4']:
        s0.append(0.)                   # He3 neutral density
        s0.append(0.)                   # He3 + ion density
        s0.append(0.)                   # He3 2+ ion density
        
    # Defining impurity related data
    for element in inp.IMPURITIES.keys():  
        # Adding impurities data
        impurity_data(element)      

        # Impurity starting condition
        if inp.IMPURITIES[element][0] == 0:
            for _ in range(ut.IMPURITIES_CHARGE_STATES[element]):
                s0.append(0.0e0)
        else:
            s0.append(inp.IMPURITIES[element][0])
            s0.append(1.0)          # Necessary for numerical stability
            for _ in range(ut.IMPURITIES_CHARGE_STATES[element]-2):
                s0.append(0.0e0)
                               
    return s0
            
# Reads fields from input report and overwrites input data          
def read_report(path: str) -> list:
    skip_rows = '-=['
    first_field = True
    fields = []
    impurities = False
    
    # Fetching all the fields data from the file
    with open(path) as f:
        for line in f:     
            if line.strip():       
                text = line.split()
                
                # Ignore everything below the last default line of the report
                if 'END' in line:
                    break
                
                # Skip unnecessary lines
                if any(char in text[0] for char in skip_rows):
                    continue
                
                # Only the first line has to be entirely store (simultation run name)
                if first_field:
                    fields.append(line.strip())
                    first_field = False
                else:
                    try:
                        fields.append(float(text[2]))
                    except Exception as e:
                        # Impurity element
                        if text[0] == 'Impurity':
                            fields.append(text[2])
                            impurities = True
                        
                        # Reactions
                        elif text[2].lower() == 'true':
                            fields.append(True)
                        elif text[2].lower() == 'false':
                            fields.append(False)
                            
                        else:
                            print(f'Exception occurred: {e}')  
    
    # Overwriting user_inputs.py data    
    # Base data
    try:
        inp.RUN_NAME = fields[0]                            # Run name   
        inp.R, inp.a = fields[1:3]                          # Reactor chamber dimentions                             
        inp.t_start, inp.t_end, inp.n_points = fields[3:6]  # Time frame
        inp.n_points = int(inp.n_points)                    # fixing n_points type to INT (otherwise it will produce an error in np.linspace)
        # Operational variables
        inp.Aux_heat_operat, inp.Aux_heat_start, inp.TimeDown = fields[6:9] # Auxiliary heating
        inp.BT, inp.IP = fields[9:11]                       # Confinement
        inp.t_ref, inp.D0_in, inp.T0_in = fields[11:14]     # Refueling   
        # Starting conditions
        base_s0 = fields[14:22]                             # Starting conditions
        # Reactions modelled
        inp.REACTIONS = {k: v for k,v in zip(inp.REACTIONS.keys(), fields[22:26])}
        
    except Exception as e:
        if e == IndexError:
            print(f'{e}.\n Some of the base parameters are missing or may not have been specified correctly')
        else:
            print(f'Exception occurred: {e}')  
    
    # Extra parameters (other than base run parameters)
    index = 26
    if len(fields) > index:
        # Impurities
        if impurities:
            if len(fields)//5 == 0:
                raise print('Impurity parameters may not be specified correctly')
            
            for i in range(len(fields[index:])//5):
                element = fields[index+i*5]
                inp.IMPURITIES[element] = fields[index+i*5+1:index+i*5+5]
                
    else:
        inp.IMPURITIES.clear()
    
    return base_s0           

# Function for adding impurities data from periodictable library and ionization potentials
def impurity_data(element: str):
    try:
        for element_object in elements:
            if element_object.symbol == element:
                ut.IMPURITIES_CHARGE_STATES[element] = element_object.number + 1
                           
                # Loading impurity ionization potentials
                with open('data/ionization_potentials','r') as f:
                    ut.POTENTIALS[element] = f.readlines()[element_object.number-1].split()
                
                # Finding the most abundant isotope for mass number
                max_abundance = 0
                chosen_isotope = None
                
                for isotope in element_object.isotopes:
                    element_isotope_object = element_object[isotope]
                    
                    if element_isotope_object.abundance is not None and element_isotope_object.abundance > max_abundance:
                        max_abundance = element_isotope_object.abundance
                        chosen_isotope = element_isotope_object
                
                if chosen_isotope:
                    ut.IMPURITIES_MASS_NUMBER[element] = chosen_isotope.mass
                else:
                    raise Exception(f'No stable isotopes found for element {element}')
                
                break
    except:
        raise Exception(f'Element {element} not recognized.')             
