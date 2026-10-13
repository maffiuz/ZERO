# %%
def main() -> dict:    
    import tools.report_manager as rm
    from tools.run_initializer import run_initializer
    from tools.numerics import numerical_solution
    from tools.plotter import plotter

    # Definition of reactor starting conditions
    s0 = run_initializer()
      
    # Numerical solution of the problem
    print('Computing the solution...')
    t, y = numerical_solution(s0)
    
    # Output
    print('Writing the report...')
    rm.startupReport()
    rm.writer(t,y)
    
    # Plot
    print('Creating the plots...')
    return plotter()

if __name__ == "__main__":
    main()
# %%
