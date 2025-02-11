function y = g1_diff(X,noise)
    tau1 = noise(1);
    E = X(3);
    y = [0;0;tau1*E;0;0];
end