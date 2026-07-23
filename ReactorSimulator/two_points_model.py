import userDefined.user_inputs as inp
from tokamak.vector import vector 
from tools.run_initializer import run_initializer
import numpy as np
import scipy.constants as cs
import matplotlib.pyplot as plt

def main():
    # Constants
    GAMMA = 7
    KAPPA0E = 2000

    # Useful data
    mass_D = 3.34e-27
    charge_D = cs.elementary_charge
    fx = 1.1

    BP = 2*1e-7*inp.IP/inp.a
    
    safety_factor_q = inp.a * inp.BT / inp.R / BP

    connection_lentgh = inp.R * cs.pi * safety_factor_q

    B_total = np.sqrt(BP**2 + inp.BT**2)

    print(BP, safety_factor_q, connection_lentgh)

    # Fetching the data
    _ = run_initializer()       # Necessary to load the dictionaries for now
    base_data = np.loadtxt("out/simulator_data",skiprows=3)

    for element in inp.IMPURITIES.keys():
        path = f'out/impurities_{element}_data'
        try:
            data = np.loadtxt(path,skiprows=5)
            base_data = np.hstack((base_data, data[:,1:]))
        except:
            raise Exception(f'Data for {element} impurity not found')

    T_target = []
    n_target = []
    Time = base_data[:,0]

    for dd in base_data:
        t = dd[0]
        values = dd[1:]

        s = vector(values)
        
        T_upstream = s.TeV * 0.01       # 1% of plasma T
        n_upstream = s.ne * 0.5
        
        P_sol = s.surfacePowerLoss(t)
        check = {}
        if P_sol < 0:
            check[t] = P_sol
        
        v_perp = np.sqrt(2*T_upstream*cs.eV/mass_D)    
        larmor_poloidal = mass_D*v_perp/charge_D/BP
        lambda_q = 2*inp.a/inp.R*larmor_poloidal
        
        q_par_upstream = P_sol/(2*cs.pi*(inp.R + inp.a)*lambda_q*(B_total/BP))
        q_par = q_par_upstream/fx
        
        T_t = pow(T_upstream**(7/2) - 7/2*q_par*connection_lentgh/KAPPA0E,2/7)
        
        T_target.append(T_t)
        n_target.append(n_upstream*T_upstream/T_t/2)
        
    print(check)

    plt.figure()
    plt.plot(Time,T_target)
    plt.grid()
    plt.xlabel('Time (s)')
    plt.ylabel('Target Temperature (eV)')

    plt.figure()
    plt.plot(Time,n_target)
    plt.grid()
    plt.xlabel('Time (s)')
    plt.ylabel('Target particle flux (m-1s-1)')

    plt.show()
