function J = jacobian_matrix(S, V, E, I, R, p)
    
    %%%%%%%%%%%%%%%%%%%%
    % Parameters
    %%%%%%%%%%%%%%%%%%%%

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the Jacobian matrix components
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    J = zeros(5,5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fill in Jacobian matrix manually (this part depends on system specifics)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dS/dt = Q - eta_s*g(I)*S - (mu + rho)*S + m*V
    J(1,1) = -eta_s * g(I,delta) - (mu + rho); % d/dS
    J(1,2) = m; % d/dV
    J(1,3) = 0; % d/dE
    J(1,4) = -eta_s * (1 + 2*delta*I) * S; % d/dI
    J(1,5) = 0; % d/dR
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dV/dt = (1-tau)*rho*S - eta_v*g(I)*V - (mu+m)*V
    J(2,1) = (1 - tau) * rho; % d/dS
    J(2,2) = -eta_v * g(I,delta) - (mu + m); % d/dV
    J(2,3) = 0; % d/dE
    J(2,4) = -eta_v * (1 + 2*delta*I) * V; % d/dI
    J(2,5) = 0; % d/dR
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dE/dt = eta_s*g(I)*S + eta_v*g(I)*V - (mu + theta)*E
    J(3,1) = eta_s * g(I,delta); % d/dS
    J(3,2) = eta_v * g(I,delta); % d/dV
    J(3,3) = -(mu + theta); % d/dE
    J(3,4) = (1 + 2*delta*I) * (eta_s * S + eta_v * V); % d/dI
    J(3,5) = 0; % d/dR
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dI/dt = theta*E - (mu + varpi)*I
    J(4,1) = 0; % d/dS
    J(4,2) = 0; % d/dV
    J(4,3) = theta; % d/dE
    J(4,4) = -(mu + varpi); % d/dI
    J(4,5) = 0; % d/dR
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dR/dt = tau*rho*S + varpi*I - mu*R
    J(5,1) = tau * rho; % d/dS
    J(5,2) = 0; % d/dV
    J(5,3) = 0; % d/dE
    J(5,4) = varpi; % d/dI
    J(5,5) = -mu; % d/dR

end