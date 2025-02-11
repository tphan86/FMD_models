%% figure 6: stochastic model - delta > delta_star, theta > theta_star, varpi = 0.12 and 0.13, small noises

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = 1000;      % recruitment rate of susceptibles
eta_s = 1E-6;  % contact rate between susceptibles and infectious
eta_v = 5E-7;  % contact rate between vaccinated and infectious
rho = 0.06;    % vaccination rate of susceptibles
tau = 0.2;    % fraction of immunized susceptibles - varied
m = 0.0056;      % immunity decline rate of vaccinated
theta = 0.05;  % progression rate from latently infected to infectious
varpi1 = 0.12;  % immunity recovery rate of infectious
varpi2 = 0.13;
delta = 1E-3; % reinfection rate of infectious - varied
mu = 0.002;    % natural mortality rate

p1 = [Q;eta_s;eta_v;rho;tau;m;theta;varpi1;delta;mu];

p2 = [Q;eta_s;eta_v;rho;tau;m;theta;varpi2;delta;mu];

tau1 = 0.05;
tau2 = 0.05;

noise = [tau1;tau2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the sde for (p1,noise) and (p2,noise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S0 = 1E6;
V0 = 1E6;
E0 = 1E3;
I0 = 1E3;
R0 = 1E5;
Xzero = [S0;V0;E0;I0;R0];
T1 = 0;
T2 = 5000;
Dt = 5*10^(-3);

% (p1,noise)
[Time1,Xrk1] = RK_stochastic_FMD(p1,noise,T1,T2,Dt,Xzero);

Xrk1_1 = [Xzero(1),Xrk1(1,:)];
Xrk1_2 = [Xzero(2),Xrk1(2,:)];
Xrk1_3 = [Xzero(3),Xrk1(3,:)];
Xrk1_4 = [Xzero(4),Xrk1(4,:)];
Xrk1_5 = [Xzero(5),Xrk1(5,:)];

% (p2,noise)
[Time2,Xrk2] = RK_stochastic_FMD(p2,noise,T1,T2,Dt,Xzero);

Xrk2_1 = [Xzero(1),Xrk2(1,:)];
Xrk2_2 = [Xzero(2),Xrk2(2,:)];
Xrk2_3 = [Xzero(3),Xrk2(3,:)];
Xrk2_4 = [Xzero(4),Xrk2(4,:)];
Xrk2_5 = [Xzero(5),Xrk2(5,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the solution paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

%%%%%%%%%%%%%%%%%

h1 = subplot(121); % varpi = 0.12
plot(Time1,Xrk1_1,'blue','LineWidth',3), hold on
plot(Time1,Xrk1_2,'magenta','LineWidth',3), hold on
plot(Time1,Xrk1_3,'cyan','LineWidth',3), hold on
plot(Time1,Xrk1_4,'red','LineWidth',3), hold on
plot(Time1,Xrk1_5,'green','LineWidth',3), 

ax = gca;
ax.FontSize = 16;
xlabel('Time (days)','FontSize',18); 
ylabel('Animals','FontSize',18);
legend('Susceptibles', 'Vaccinated', 'Exposed',...
       'Infectious', 'Recovered', ...
       'FontSize', 14,'Location', 'Best');
ax.TitleHorizontalAlignment = 'left';
title('$\textbf{A.}$ $\lambda>0$, $\varpi = 0.12$',...
      'Interpreter','latex','FontSize',20)
ytickformat('%.2f')
xlim([0 3000])
set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%

h2 = subplot(122); % varpi = 0.14
plot(Time2,Xrk2_1,'blue','LineWidth',3), hold on
plot(Time2,Xrk2_2,'magenta','LineWidth',3), hold on
plot(Time2,Xrk2_3,'cyan','LineWidth',3), hold on
plot(Time2,Xrk2_4,'red','LineWidth',3), hold on
plot(Time2,Xrk2_5,'green','LineWidth',3), 

ax = gca;
ax.FontSize = 16;
xlabel('Time (days)','FontSize',18); 
ylabel('Animals','FontSize',18);
legend('Susceptible', 'Vaccinated', 'Exposed',...
       'Infectious', 'Recovered', ...
       'FontSize', 14,'Location', 'Best');
ax.TitleHorizontalAlignment = 'left';
title('$\textbf{B.}$ $\lambda>0$, $\varpi = 0.13$',...
      'Interpreter','latex','FontSize',20)
ytickformat('%.2f')
xlim([0 3000])
set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print out the plot as a pdf file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'Units','normalized')

set(h2,'position',[.58 .13 .40 .75])
set(h1,'position',[.06 .13 .40 .75])
set(gcf, 'PaperSize', [15 8], 'PaperPosition', [0 0 15 8])
print('figure6','-dpdf')