import numpy as np
from scipy.integrate import solve_ivp

import userDefined.user_inputs as inp
import tools.utilities as ut
from tools.auxiliary import Evolution

def numerical_solution(s0: list) -> tuple[list, list]:
    # Sorting impurities for their activation times
    ut.IMPURITIES_SOLUTION = {k: v for k, v in sorted(ut.IMPURITIES_SOLUTION.items(), key=lambda item: item[1][0])}
    
    # Definition of the time evolution function of vector s
    F = lambda t, s: Evolution(t,s)  
    
    # Building the solution logic control
    solution = []   # [t_span, impurities]
    imp = [] 
    
    # If there are no impurities active at t=0
    if not any(active for _,active in ut.IMPURITIES_SOLUTION.values()):
        # Check for any impurities
        if ut.IMPURITIES_SOLUTION:
            t_end,_ = next(iter(ut.IMPURITIES_SOLUTION.values()))
            t_end = min(t_end,inp.t_end)
        else:
            t_end = inp.t_end     

    else:
        # Looking for all active impurities at t=0
        for element, (activation_time, active) in ut.IMPURITIES_SOLUTION.items():
            if activation_time == 0:
                imp.append(element)
            else:
                t_end = activation_time
                
    # Computing first step solution
    t_span, t_glob, y_glob, s0_next = solver(F, inp.t_start, t_end, s0)            
    solution.append([t_span, imp.__repr__()]) 
   
    # Defining impurity segments
    t_start = t_end
    for element, (activation_time, active) in ut.IMPURITIES_SOLUTION.items():
        if activation_time == 0:
            continue
        
        if activation_time >= inp.t_end:
            break
        
        if activation_time > t_start:
            t_end = activation_time
            t_span, t, y, s0_next = solver(F, t_start, t_end, s0_next) 
            solution.append([t_span, imp.__repr__()]) 
            
            t_glob = np.concatenate([t_glob, t[1:]])
            y_glob = np.concatenate([y_glob, y[:,1:]],axis=1)  
        
            t_start = t_end
        
        # For impurities not present at t=0, expand the s0 and add zeros in the previous solution
        imp.append(element)
        ut.IMPURITIES_SOLUTION[element][1] = True
        s0_next = np.concatenate([s0_next, np.zeros(ut.IMPURITIES_CHARGE_STATES[element])])
        y_glob = np.vstack([y_glob, np.zeros((ut.IMPURITIES_CHARGE_STATES[element],y_glob.shape[1]))])
    
    # Last solution segment 
    if t_start < inp.t_end:
        t_span, t, y, s0_next = solver(F, t_start, inp.t_end, s0_next) 
        solution.append([t_span, imp.__repr__()]) 
        
        t_glob = np.concatenate([t_glob, t[1:]])
        y_glob = np.concatenate([y_glob, y[:,1:]],axis=1)  
    
    # Small log info
    print('\nSOLUTION STRUCTURE')
    for row in solution:
        print(f't_span: {row[0][0]} s to {row[0][1]} s\t Impurities: {row[1]}')  
        
    print('')  
    
    return t_glob, y_glob  


def solver(F: function, t_start: float, t_end: float, s0: list) -> tuple[list, list, list, list]:
    t_span = [t_start, t_end]
    n_points = int(inp.n_points*(t_end-t_start)/(inp.t_end-inp.t_start))
    t_eval = np.linspace(t_start, t_end, n_points)
        
    # Computing time evolution with Implicit Runge-Kutta method of the Radau IIA family of order 5 
    sol = solve_ivp(F, t_span, s0, method = 'Radau', t_eval=t_eval, dense_output=True)
    if not sol.success:
        raise Exception(f'The solver did not reach convergence. Exception raised: {sol.message}')
    print('--- Solved timespan')
    t = sol.t
    y = sol.y
    s0_next = sol.sol(t_end)    
    
    return t_span, t, y, s0_next
