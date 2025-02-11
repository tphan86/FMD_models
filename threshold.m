function y = threshold(p,noise)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    tau1 = noise(1); % noise intensity for exposed
    tau2 = noise(2); % noise intensity for infectious

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % combined parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S0 = Q * (mu + m) / (mu^2 + (m + rho) * mu + m * tau * rho);
    V0 = Q * (1 - tau) * rho / (mu^2 + (m + rho) * mu + m * tau * rho);

    theta_bar = theta*(eta_s*S0 + eta_v*V0);

    w = 4*sqrt(theta_bar) / (tau1^2 + tau2^2);

    Theta = 2*(varpi + tau2^2 - theta) / (tau1^2 + tau2^2) - 1;

    R_Theta = besselk(Theta + 1, w) / besselk(Theta, w);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % threshold lambda
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    y = sqrt(theta_bar)*R_Theta - mu - varpi - tau2^2/2;
    
end