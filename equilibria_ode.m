function [DFE, EE] = equilibria_ode(p)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unpack parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Q = p(1);      % recruitment rate of susceptibles
    eta_s = p(2);  % contact rate between susceptibles and infectious
    eta_v = p(3);  % contact rate between vaccinated and infectious
    rho = p(4);    % vaccination rate of susceptibles
    tau = p(5);    % fraction of immunized susceptibles
    m = p(6);      % immunity decline rate of vaccinated
    theta = p(7);  % progression rate from latently infected to infectious
    varpi = p(8);  % immunity recovery rate of infectious
    delta = p(9); % reinfection rate of infectious
    mu = p(10);    % natural mortality rate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Disease Free Equilibrium (DFE)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S0 = Q * (mu + m) / (mu^2 + (m + rho) * mu + m * tau * rho);
    V0 = Q * (1 - tau) * rho / (mu^2 + (m + rho) * mu + m * tau * rho);
    R0 = Q * rho * tau * (mu + m) / (mu^3 + (m + rho) * mu^2 + m * tau * rho * mu);

    DFE = [S0, V0, 0, 0, R0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coefficients for the quartic equation in I (for Endemic Equilibrium)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A = (mu + theta) * (mu + varpi) / theta;
    w1 = A * eta_v * eta_s;
    w2 = A * (eta_s * (mu + m) + eta_v * (mu + rho));
    w3 = A * (mu^2 + (m + rho) * mu + m * tau * rho);
    w4 = Q * eta_v * eta_s;
    w5 = Q * eta_s * (mu + m);
    w6 = Q * eta_v * (1 - tau) * rho;
    R_0 = (eta_s*S0 + eta_v*V0) / A; % Reproduction number
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the quartic equation in I
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    syms I;
    f = w1 * delta^2 * I^4 + (2 * delta * w1 - delta^2 * w4) * I^3 + ...
        (w1 + delta * w2 - 2 * delta * w4) * I^2 + ...
        (w2 - w4 - delta * w5 - delta * w6) * I + w3 * (1 - R_0) == 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve the quartic equation for I
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    I_solutions = double(vpasolve(f, I));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filter to keep only positive, real-valued solutions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    I_solutions = I_solutions(imag(I_solutions) == 0 & I_solutions > 0); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the endemic equilibrium for each solution of I
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    EE = [];

    for I = I_solutions'

        S = Q * (eta_v * delta * I^2 + eta_v * I + mu + m) / ...
            ((eta_s * delta * I^2 + eta_s * I + mu + rho) *...
            (eta_v * delta * I^2 + eta_v * I + mu + m) - m * (1 - tau) * rho);

        V = Q * (1 - tau) * rho / ...
            ((eta_s * delta * I^2 + eta_s * I + mu + rho) *...
            (eta_v * delta * I^2 + eta_v * I + mu + m) - m * (1 - tau) * rho);

        E = (mu + varpi) * I / theta;

        R = tau * rho * S / mu + varpi * I / mu;

        EE = [EE; S, V, E, I, R];
    end
    
    %%%%%%%%%%%%%%%%%%%%
    % Display results
    %%%%%%%%%%%%%%%%%%%%

    fprintf('The reproduction number: %.4f\n', R_0)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Disease Free Equilibrium:\n');
    fprintf('S0 = %.4f, V0 = %.4f, E0 = 0, I0 = 0, R0 = %.4f\n', S0, V0, R0);
    
    % Jacobian at DFE
    J_DFE = jacobian_matrix(S0, V0, 0, 0, R0, p);
    eigen_DFE = eig(J_DFE);
    
    % Check stability of DFE
    if all(real(eigen_DFE) < 0)
        fprintf('The DFE is stable (all eigenvalues have negative real parts).\n');
    else
        fprintf('The DFE is unstable (at least one eigenvalue has positive real part).\n');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nEndemic Equilibria (EE):\n');

    if size(EE,1) == 0

        fprintf('None\n')

    else
    
        for i = 1:size(EE, 1)
            
            fprintf('EE%d: S = %.4f, V = %.4f, E = %.4f, I = %.4f, R = %.4f\n', ...
                    i, EE(i, 1), EE(i, 2), EE(i, 3), EE(i, 4), EE(i, 5));

            % Jacobian at EE
            J_EE = jacobian_matrix(EE(i, 1), EE(i, 2), EE(i, 3), EE(i, 4), EE(i, 5), p);
            eigen_EE = eig(J_EE);
            
            % Check stability of EE
            if all(real(eigen_EE) < 0)
                fprintf('EE%d is stable (all eigenvalues have negative real parts).\n', i);
            else
                fprintf('EE%d is unstable (at least one eigenvalue has positive real part).\n', i);
            end
    
        end

    end

end