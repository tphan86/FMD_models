function y = f_drift(X,p)

    %%%%%%%%%%%%%%%%%%%%%
    % parameters
    %%%%%%%%%%%%%%%%%%%%%

    Q = p(1);      % recruitment rate of susceptibles
    eta_s = p(2);  % contact rate between susceptibles and infectious
    eta_v = p(3);  % contact rate between vaccinated and infectious
    rho = p(4);    % vaccination rate of susceptibles
    tau = p(5);    % fraction of immunized susceptibles
    m = p(6);      % immunity decline rate of vaccinated
    theta = p(7);  % progression rate from latently infected to infectious
    varpi = p(8);  % immunity recovery rate of infectious
    delta = p(9);  % reinfection rate of infectious - convex parameter
    mu = p(10);    % natural mortality rate
    
    %%%%%%%%%%%%%%%%%%%%%
    % variables
    %%%%%%%%%%%%%%%%%%%%%
    
    S = X(1);  % susceptible population
    V = X(2);  % vaccinated population
    E = X(3);  % exposed (latently infected) population
    I = X(4);  % infectious population
    R = X(5);  % recovered population

    %%%%%%%%%%%%%%%%%%%%%
    % Equations
    %%%%%%%%%%%%%%%%%%%%%
    
    dS = Q - eta_s*g(I,delta).*S - (mu + rho)*S + m*V;
    dV = (1 - tau)*rho*S - eta_v*g(I,delta).*V - (mu + m)*V; 
    dE = eta_s*g(I,delta).*S + eta_v*g(I,delta).*V - (mu + theta)*E; 
    dI = theta*E - (mu + varpi)*I; 
    dR = tau*rho*S + varpi*I - mu*R;

    %%%%%%%%%%%%%%%%%%%%%%
    % Output
    %%%%%%%%%%%%%%%%%%%%%%

    y = [dS;dV;dE;dI;dR]; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end