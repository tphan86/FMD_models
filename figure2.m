%% figure 2: delta > delta_star, theta > theta_star, small noises

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
varpi = 0.1;  % immunity recovery rate of infectious
delta = 1E-3; % reinfection rate of infectious - varied
mu = 0.002;    % natural mortality rate

p = [Q;eta_s;eta_v;rho;tau;m;theta;varpi;delta;mu];

tau1 = 0.05;
tau2 = 0.05;

noise = [tau1;tau2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the ode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DFE, EE] = equilibria_ode(p);

check_stability_sde(p,noise);

func = @f_drift;

S0 = round(EE(1),2) + 1E4;
V0 = round(EE(2),2) + 1E4;
E0 = round(EE(3),2) + 1E4;
I0 = round(EE(4),2) + 1E4;
R0 = round(EE(5),2) + 1E4;

initial_value = [S0;V0;E0;I0;R0];

tspan = [0,5000];

[t1,y1] = ode45(@(t,y) func(y,p),tspan,initial_value);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the sde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T1 = 0;
T2 = 5000;
Dt = 5*10^(-3);
Xzero = initial_value;
rng(250)
[t2,Xrk2] = RK_stochastic_FMD(p,noise,T1,T2,Dt,Xzero);

y2_1 = [Xzero(1),Xrk2(1,:)];
y2_2 = [Xzero(2),Xrk2(2,:)];
y2_3 = [Xzero(3),Xrk2(3,:)];
y2_4 = [Xzero(4),Xrk2(4,:)];
y2_5 = [Xzero(5),Xrk2(5,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the solution paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

%%%%%%%%%%%%%%%%%

h1 = subplot(221); % ODE solutions

plot(t1,y1(:,1),'blue','LineWidth',3), hold on
plot(t1,y1(:,2),'magenta','LineWidth',3), hold on
plot(t1,y1(:,3),'green','LineWidth',3),

ax = gca;
ax.FontSize = 16;

ylabel('Animals','FontSize',18);
legend('Susceptible', 'Vaccinated', 'Removed',...
       'FontSize', 14,'Location', 'Best');
ax.TitleHorizontalAlignment = 'left';
title('$\textbf{A.}$ $\delta>\delta^*$, $\theta>\theta^*$, $\tau_1=\tau_2=0$',...
      'Interpreter','latex','FontSize',20)

axis([0 3000 1E3 1E5])
set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%

h2 = subplot(222); % SDE solutions

plot(t2,y2_1,'blue','LineWidth',3), hold on
plot(t2,y2_2,'magenta','LineWidth',3), hold on
plot(t2,y2_5,'green','LineWidth',3),

ax = gca;
ax.FontSize = 16;

ylabel('Animals','FontSize',18);
legend('Susceptible', 'Vaccinated', 'Removed', ...
       'FontSize', 14,'Location', 'Best');
ax.TitleHorizontalAlignment = 'left';
title('$\textbf{B.}$ $\delta>\delta^*$, $\theta>\theta^*$, $\tau_1 = \tau_2 = 0.05$',...
      'Interpreter','latex','FontSize',20)

xlim([0 3000])
set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%%%%%

h3 = subplot(223); % ODE solutions

plot(t1,y1(:,3),'cyan','LineWidth',3), hold on
plot(t1,y1(:,4),'red','LineWidth',3), 

ax = gca;
ax.FontSize = 16;
xlabel('Time (days)','FontSize',18); 
ylabel('Animals','FontSize',18);
legend('Exposed', 'Infectious', ...
       'FontSize', 14,'Location', 'Best');
ytickformat('%.2f')
axis([0 3000 1E3 1E5])
set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%

h4 = subplot(224); % SDE solution

plot(t2,y2_3,'cyan','LineWidth',3), hold on
plot(t2,y2_4,'red','LineWidth',3),

ax = gca;
ax.FontSize = 16;
xlabel('Time (days)','FontSize',18); 
ylabel('Animals','FontSize',18);
legend('Exposed', 'Infectious', ...
       'FontSize', 14,'Location', 'Best');
ytickformat('%.2f')
xlim([0 3000])
set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print out the plot as a pdf file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'Units','normalized')
set(h4,'position',[.58 .1 .40 .38])
set(h3,'position',[.06 .1 .40 .38])
set(h2,'position',[.58 .54 .40 .38])
set(h1,'position',[.06 .54 .40 .38])
set(gcf, 'PaperSize', [15 10], 'PaperPosition', [0 0 15 10])
print('figure2','-dpdf')