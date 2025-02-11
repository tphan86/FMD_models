function [Time,Xrk] = RK_stochastic_FMD(p,noise,T1,T2,Dt,Xzero)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preparation: step size, number of iterations for RK method,
    % number of iterations for EM method for Brownian Increments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % rng(250) % random number generator

    T = T2 - T1; % amount of time for integration
    n = T/Dt;    % number of iterations for RK method
                 % note: Dt = time step for RK method
    dt = (Dt)^2; % time step for EM method
    N = Dt/dt;   % number of iterations for EM method
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pre-allocate matrices for RK methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Time = zeros(1,n+1); % store time points
    Time(1) = T1;
    Xrk = zeros(5,n);    % store solution values at time points
    Xtilde = zeros(5,5); % store supporting values
    Xfinal = Xzero;      % store solution values after each iteration
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop through the number of iterations for RK method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for k = 1:n

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EM method for dW21, dW2, dW12, and dW1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Y1 = 0; Y2 = 0; 
        Z1 = 0; Z2 = 0;

        for m = 1:N
            dW1 = sqrt(dt)*randn;
            dW2 = sqrt(dt)*randn;
            Y1 = Y1 + Y2*dW1;
            Y2 = Y2 + dW2;
            Z1 = Z1 + Z2*dW2;
            Z2 = Z2 + dW1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Runge-Kutta method for the stochastic FMD model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        X = Xfinal;

        Xtilde(:,1) = X + f_drift(X,p)*Dt;

        Xtilde(:,2) = Xtilde(:,1) + g1_diff(X,noise)*0.5*(Z2^2-Dt)/sqrt(Dt)...
                                  + g2_diff(X,noise)*Y1/sqrt(Dt);

        Xtilde(:,3) = Xtilde(:,1) + g1_diff(X,noise)*Z1/sqrt(Dt)...
                                  + g2_diff(X,noise)*0.5*(Y2^2-Dt)/sqrt(Dt);

        Xtilde(:,4) = Xtilde(:,1) - g1_diff(X,noise)*0.5*(Z2^2-Dt)/sqrt(Dt)...
                                  - g2_diff(X,noise)*Y1/sqrt(Dt);

        Xtilde(:,5) = Xtilde(:,1) - g1_diff(X,noise)*Z1/sqrt(Dt)...
                                  - g2_diff(X,noise)*0.5*(Y2^2-Dt)/sqrt(Dt);

        Xfinal = X + 0.5*(f_drift(X,p)+f_drift(Xtilde(:,1),p))*Dt...
                   + Z2*g1_diff(X,noise)... 
                   + 0.5*sqrt(Dt)*(g1_diff(Xtilde(:,2),noise)-g1_diff(Xtilde(:,4),noise))...
                   + Y2*g2_diff(X,noise)... 
                   + 0.5*sqrt(Dt)*(g2_diff(Xtilde(:,3),noise)-g2_diff(Xtilde(:,5),noise));

        Xrk(:,k) = Xfinal;

        Time(k+1) = k*Dt;

    end
end