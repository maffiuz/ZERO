# %%
def main():
    import numpy as np
    from scipy.integrate import solve_ivp
    
    import userDefined.user_inputs as inp
    import tools.report_manager as rm
    from tools.run_initializer import run_initializer
    from tools.auxiliary import Evolution
    from tools.plotter import plotter

    # Definition of reactor and starting conditions
    s0 = run_initializer()
    
    # Definition of the time evolution function of vector s
    F = lambda t, s: Evolution(t,s)  
    
    # Defining the time projection of the solution
    t_eval = np.linspace(inp.t_start, inp.t_end, inp.n_points)
        
    # Computing time evolution with Implicit Runge-Kutta method of the Radau IIA family of order 5
    sol = solve_ivp(F, [inp.t_start, inp.t_end], s0, method = 'Radau', t_eval=t_eval)
    
    if not sol.success:
        raise Exception(f'The solver did not reach convergence. Exception raised: {sol.message}')
    
    # Output
    rm.startupReport()
    rm.writer(sol.t,sol.y)
    
    # Plot
    plotter()

if __name__ == "__main__":
    main()
# %%
