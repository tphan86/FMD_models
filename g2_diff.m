function y = g2_diff(X,noise)
    tau2 = noise(2);
    I = X(4);
    y = [0;0;0;tau2*I;0];
end