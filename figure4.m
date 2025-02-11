%% figure 4: deterministic model - delta > delta_star, theta > theta_star, varpi = 0.12 and 0.13

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the ode with p1 and p2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DFE1, EE1] = equilibria_ode(p1);
[DFE2, EE2] = equilibria_ode(p2);
func = @f_drift;

S0_1 = round(EE1(1),2) + 1E2;
V0_1 = round(EE1(2),2) + 1E2;
E0_1 = round(EE1(3),2) + 1E0;
I0_1 = round(EE1(4),2) + 1E0;
R0_1 = round(EE1(5),2) + 1E2;
initial_value_1 = [S0_1;V0_1;E0_1;I0_1;R0_1];

S0_2 = round(EE2(1),2) + 1E2;
V0_2 = round(EE2(2),2) + 1E2;
E0_2 = round(EE2(3),2) + 1E0;
I0_2 = round(EE2(4),2) + 1E0;
R0_2 = round(EE2(5),2) + 1E2;
initial_value_2 = [S0_2;V0_2;E0_2;I0_2;R0_2];

tspan = [0,5000];

[t1,y1] = ode45(@(t,y) func(y,p1),tspan,initial_value_1);
[t2,y2] = ode45(@(t,y) func(y,p2),tspan,initial_value_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the solution paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

%%%%%%%%%%%%%%%%%

h1 = subplot(121); % varpi1 = 0.12

plot(t1,y1(:,1),'blue','LineWidth',3), hold on
plot(t1,y1(:,2),'magenta','LineWidth',3), hold on
plot(t1,y1(:,3),'cyan','LineWidth',3), hold on
plot(t1,y1(:,4),'red','LineWidth',3), hold on
plot(t1,y1(:,5),'green','LineWidth',3),

ax = gca;
ax.FontSize = 16;
xlabel('Time (days)','FontSize',18); 
ylabel('Animals','FontSize',18);
legend('Susceptible', 'Vaccinated', 'Exposed',...
       'Infectious', 'Recovered', ...
       'FontSize', 14,'Location', 'Best');
ax.TitleHorizontalAlignment = 'left';
title('$\textbf{A.}$ $\delta>\delta^*$, $\theta>\theta^*$, $\varpi=0.12$',...
      'Interpreter','latex','FontSize',20)
ytickformat('%.2f')
ylim([1E3 1E6])
set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%

h2 = subplot(122); % varpi2 = 0.14

plot(t2,y2(:,1),'blue','LineWidth',3), hold on
plot(t2,y2(:,2),'magenta','LineWidth',3), hold on
plot(t2,y2(:,3),'cyan','LineWidth',3), hold on
plot(t2,y2(:,4),'red','LineWidth',3), hold on
plot(t2,y2(:,5),'green','LineWidth',3),

ax = gca;
ax.FontSize = 16;
xlabel('Time (days)','FontSize',18); 
ylabel('Animals','FontSize',18);
legend('Susceptible', 'Vaccinated', 'Exposed',...
       'Infectious', 'Recovered', ...
       'FontSize', 14,'Location', 'Best');
ax.TitleHorizontalAlignment = 'left';
title('$\textbf{B.}$ $\delta>\delta^*$, $\theta>\theta^*$, $\varpi=0.13$',...
      'Interpreter','latex','FontSize',20)
ytickformat('%.2f')
set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print out the plot as a pdf file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'Units','normalized')

set(h2,'position',[.58 .13 .40 .75])
set(h1,'position',[.06 .13 .40 .75])
set(gcf, 'PaperSize', [15 8], 'PaperPosition', [0 0 15 8])
print('figure4','-dpdf')