function lambda = lyapunov_exponents(p,noise,varpi)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    p(8) = varpi; % Note that at the beginning p(8) = 0

    S0 = 1E6;
    V0 = 1E6;
    E0 = 1E3;
    I0 = 1E3;
    R0 = 1E5;
    Xzero = [S0;V0;E0;I0;R0];
    T1 = 0;
    T2 = 5000;
    T = T2 - T1;
    Dt = 5*10^(-3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute SDE solutions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Time, Xrk] = RK_stochastic_FMD(p,noise,T1,T2,Dt,Xzero);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute Lyapunov exponents
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lambda = [log(Xrk(1,end))/T 
              log(Xrk(2,end))/T 
              log(Xrk(3,end))/T
              log(Xrk(4,end))/T
              log(Xrk(5,end))/T];

end