import numpy as np
import scipy.constants as cs
from scipy.interpolate import RectBivariateSpline as coeffInterp
from numpy.linalg import norm as norm

def get_next_field(istart: int, linestring: str) -> tuple[str,int]:
    # Function used to get single data fields
    no_field = [' ','\n',':'] 
    
    thischar = linestring[istart]
    while thischar in no_field:
        istart += 1 
        thischar = linestring[istart]
    # Get the start
    
    iend = istart
    thischar = linestring[iend]
    while not thischar in no_field:
        iend += 1 
        thischar = linestring[iend]
    # Get the end
    
    return linestring[istart:iend],iend

def get_next_array(istart,linestring,number_of_data):
    # Function used to get an array of length number_of_data of elements using the get_next_field function
    
    tmp = []
    for ntok in range(number_of_data):
        strfield,istart = get_next_field(istart,linestring)
        tmp.append(float(strfield))
    istart = next_line(istart,linestring)   
    return tmp,istart

def next_line(istart,linestring):
    # Function to move the reader variable to the next line
    
    while not linestring[istart]=='\n':
        istart += 1
    return istart+1

def get_next_table(istart,linestring,NROWS,NCOLS):
    # Function used to get an array of lenght NROWS of elements using the get_next_array function
    # The number of rows NROWS it's used to iterate the get_next_array function
    
    output = []
    for ITEMP in range(NROWS):
        tmp,istart = get_next_array(istart,linestring,NCOLS)
        output.append(tmp)
        # Read a single table row

    return output,istart                

def convert_to_MKS_eV(rate):
    # Convert units to MKS_eV  
    L10cm3_to_L10m3 = 6
    
    for ID in range(rate['IDMAX']):
        rate['DDENSD'][ID] += L10cm3_to_L10m3
        
    for Z1 in range(rate['IZ1MIN'],rate['IZ1MAX']+1):
        this_charge_state = 'Z='+str(Z1)
        for ID in range(rate['IDMAX']):
            for IT in range(rate['ITMAX']):
                rate['DRCOFD'][this_charge_state][IT][ID] -= L10cm3_to_L10m3
       
    return rate
    
def reads_acd(file_path):
    """Function created to read the files acd (& scd) downloaded from the OPEN-ADAS project.
    We have preferred to create our own function for the reading of the data rather than using
    a premade one. The files contain the values of the rate coefficients of ionisation/recombination
    of the elements considered.

    Args:
        file_path (string): the file path of the acd (or scd) .dat
        
    Returns:
        dictionary: the function builds a dictionary with the following keys
        - IZMAX (int): the number of different ionization states
        - IDMAX (int): the number of density values (log10) [m-3] the table is built upon
        - ITMAX (int): the number of absolute temperature values (log10) [K] the table is built upon
        - IZ1MIN (int): the minimum ionization state
        - IZ1MAX (int): the maximum ionization state
        - DDENSD (array (float)): the log10 values of density considered
        - DTEVD (array (float)): the log10 values of absolute temperature considered
        - DRCOFD (dictionary): contains values for every charge state named "Z=n". Those values are stored in 
            tables Density-Temperature
        - EVALUATE (class): provides an interpolation class to properly output a value, given any condition (whithin bounds)
    """
    
    with open(file_path, 'r') as file:
        file_content = file.read()
    # Reads the datafile
    
    output = {}
    istart = 0
    strfield,istart = get_next_field(istart,file_content)
    output['IZMAX'] = int(strfield)
    strfield,istart = get_next_field(istart,file_content)
    output['IDMAX'] = int(strfield)
    strfield,istart = get_next_field(istart,file_content)
    output['ITMAX'] = int(strfield)
    strfield,istart = get_next_field(istart,file_content)
    output['IZ1MIN'] = int(strfield)
    strfield,istart = get_next_field(istart,file_content)
    output['IZ1MAX'] = int(strfield)
    istart = next_line(istart,file_content)
    istart = next_line(istart,file_content)
    # Reads the dimensions
    
    tmp,istart = get_next_array(istart,file_content,output['IDMAX'])
    output['DDENSD'] = tmp
    tmp,istart = get_next_array(istart,file_content,output['ITMAX'])
    output['DTEVD'] = tmp
    # Reads the density and temperature arrays
    
    charge_state = {}
    for Z1 in range(output['IZ1MIN'],output['IZ1MAX']+1):
        istart = next_line(istart,file_content)
        this_charge_state = 'Z='+str(Z1)
        # charge_state[this_charge_state],istart = get_next_array(istart,file_content,output['IDMAX'])
# NEED TO TEST THAT THE FOLLOWING SUBSTITUTION IS HARMLESS        
        charge_state[this_charge_state],istart = get_next_table(istart,file_content,output['ITMAX'],output['IDMAX'])
        # charge_state[this_charge_state] = []
        # for ITEMP in range(output['ITMAX']):
        #     tmp,istart = get_next_array(istart,file_content,output['IDMAX'])
        #     charge_state[this_charge_state].append(tmp)
            # Read coefficients for charge state Z1                
    output['DRCOFD'] = charge_state
    
    output = convert_to_MKS_eV(output)
    
    evaluate = {}
    for Z1 in range(output['IZ1MIN'],output['IZ1MAX']+1):
        this_charge_state = 'Z='+str(Z1)
        # evaluate[this_charge_state] = []
        evaluate[this_charge_state] = coeffInterp(output['DTEVD'],output['DDENSD'],output['DRCOFD'][this_charge_state])
    output['EVALUATE'] = evaluate
    # Also, provide a proper interpolation function
        
    return output

def reads_individual_line(istart,file_content):
    """Function created to read the individual lines of pec files downloaded from the OPEN-ADAS project.

    Args:
        istart (int): the current reading starting point
        file_content (array): the content of the file
        
    Returns:
        dictionary: the function builds a dictionary with the following keys
        - WAVE_LENGTH (str): the color of the emission
        - NDENS (int): the number of density values (log10) considered 
        - NTEMP (int): the number of absolute temperature values (log10) considered
        - TYPE (str)
        - DENS (array (float)): the log10 values of density considered
        - TEMP (array (float)): the log10 values of absolute temperature considered
        - PEC (array (array (float))): the photon emissivity coefficients
    """  
    
    output = {}

    strfield,istart = get_next_field(istart,file_content)    
    output['WAVE_LENGTH'] = strfield

    strfield,istart = get_next_field(istart,file_content)
    # Sometimes the unit A is separated by a space
    try:
        int(strfield)
    except:
        strfield,istart = get_next_field(istart,file_content)
        
    output['NDENS'] = int(strfield)

    strfield,istart = get_next_field(istart,file_content)
    output['NTEMP'] = int(strfield)
    
    # _,istart = get_next_field(istart,file_content)
    # _,istart = get_next_field(istart,file_content)
    # _,istart = get_next_field(istart,file_content)
    # _,istart = get_next_field(istart,file_content)
    # _,istart = get_next_field(istart,file_content)
    # strfield,istart = get_next_field(istart,file_content)
    # output['TYPE'] = strfield
    
    istart = next_line(istart,file_content)
    
    tmp,istart = get_next_array(istart,file_content,output['NDENS'])
    output['DENS'] = tmp
    
    tmp,istart = get_next_array(istart,file_content,output['NTEMP'])
    output['TEMP'] = tmp
    
    tmp,istart = get_next_table(istart,file_content,output['NDENS'],output['NTEMP'])
    output['PEC'] = tmp
    
    return output,istart

def reads_pec(file_path: str) -> dict:  
    """Function created to read the files pec downloaded from the OPEN-ADAS project.
    We have preferred to create our own function for the reading of the data rather than using
    a premade one. This files contain the values of photon emissivity coefficients per wawelenght
    of the elements considered.

    Args:
        file_path (string): the file path of the pec .dat
        
    Returns:
        dictionary: the function builds a dictionary with the following keys
        - NUMBER_OF_LINES (int): the number of different wavelenghts considered
        - SPECIES (str): the atomic symbol of the species
        - CHARGE_STATE (str): the charge state of the species
        - INDIVIDUAL_LINES (array (dictionary)): array of all the individual lines dictionaries built
            using the reads_individual_line function
        - REC (class): provides an interpolation class to properly output a value, given any condition (whithin bounds)
    """  

    # Reads the datafile
    with open(file_path, 'r') as file:
        file_content = file.read()
    
    # Initializating reader
    output = {}
    istart = 0
    
    # Reads the preambule
    strfield,istart = get_next_field(istart,file_content)
    output['NUMBER_OF_LINES'] = int(strfield)
    
    strfield,istart = get_next_field(istart,file_content)
    # Adjusting formatting for atomic element
    strfield = strfield.strip('/+').lower().capitalize()
    output['SPECIES'] = strfield
    
    strfield,istart = get_next_field(istart,file_content)
    # ADAS formatting isn't consistent, therefore it's necessary to check if there is a space
    # so H + counts as 2 fields
    # or there isn't, like He+
    try:
        int(strfield)
    except:
        strfield,istart = get_next_field(istart,file_content)
    output['CHARGE_STATE'] = int(strfield)

    istart = next_line(istart,file_content)
    
    # In some PEC files there is an extra preambule that needs to be skipped
    while True:
        strfield,istart_new = get_next_field(istart,file_content)
        try:
            if float(strfield.strip('A'))%1 != 0:
                # Exit if a floating point number is found (first wavelength)
                # Removing A is crucial since in some datasets the unit is attached to the number
                break
        except:
            pass
        # Updating istart
        istart = istart_new
          
    # Reads all the spectral lines
    output['INDIVIDUAL_LINES'] = []
    for _ in range(output['NUMBER_OF_LINES']):
        newline,istart = reads_individual_line(istart,file_content)
        output['INDIVIDUAL_LINES'].append(newline)
    
    
    reference_line = output['INDIVIDUAL_LINES'][0]
    length = float(reference_line['WAVE_LENGTH'].replace('A','0'))*1.e-10
    energy = cs.h*cs.c/length
    REC = np.array(reference_line['PEC'])*energy
    # print('Loaded line at '+reference_line['WAVE_LENGTH']+' of type '+reference_line['TYPE'])
    # Evaluate REC for the reference line
    
    for single_line in output['INDIVIDUAL_LINES'][1:]:
        # if norm(np.subtract(reference_line['DENS'],single_line['DENS'])) > norm(reference_line['DENS']):
        #     print('Density error importing line '+single_line['WAVE_LENGTH']+': '+single_line['TYPE'])
        # if norm(np.subtract(reference_line['TEMP'],single_line['TEMP'])) > norm(reference_line['TEMP']):
        #     print('Temperature error importing line '+single_line['WAVE_LENGTH']+': '+single_line['TYPE'])
        length = float(reference_line['WAVE_LENGTH'].replace('A','0'))*1.e-10
        energy = cs.h*cs.c/length
        REC += np.array(single_line['PEC'])*energy
        # print('Loaded line at '+single_line['WAVE_LENGTH']+' of type '+single_line['TYPE'])
    
    cm3_to_m3 = 1.e6
    REC = np.log10(np.transpose(REC/cm3_to_m3))
    L10Te = np.log10(single_line['TEMP'])
    L10ne = np.log10(single_line['DENS'])
    # Convert units from cm3 s-1 to m3 s-1
    # Transpose to make force the first index to correspond to temperatures
    # Make sure REC is given in logarithmic scale
      
    output['REC'] = coeffInterp(L10Te,L10ne,REC)
    return output
