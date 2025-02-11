function check_stability_sde(p,noise)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the threshold
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    lambda = threshold(p,noise);
    
    %%%%%%%%%%%%%%%%%
    % Check the stability of PI_0 and PI_star
    if lambda < 0
        fprintf('threshold lambda = %.4f\n', lambda)
        fprintf('The solution of the stochastic FMD system converges a.s. to DFE\n')
    end
    if lambda > 0
        fprintf('threshold lambda = %.4f\n', lambda)
        fprintf(['The stochastic FMD system is stochastically persistent' ...
                 'in the sense that its solution converges weakly to a unique ' ...
                 'invariant probability PI* (the endemic equilibrium state)'])
    end
    %%%%%%%%%%%%%%%%
end