%% figure 5: Lyapunov exponents when lambda > 0 and varpi changes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = 1000;      % recruitment rate of susceptibles
eta_s = 1E-6;  % contact rate between susceptibles and infectious
eta_v = 5E-7;  % contact rate between vaccinated and infectious
rho = 0.06;    % vaccination rate of susceptibles
tau = 0.2;    % fraction of immunized susceptibles - varied
m = 0.0056;      % immunity decline rate of vaccinated
theta = 0.05;  % progression rate from latently infected to infectious
varpi = 0;  % immunity recovery rate of infectious
delta = 1E-3; % reinfection rate of infectious - varied
mu = 0.002;    % natural mortality rate

p = [Q;eta_s;eta_v;rho;tau;m;theta;varpi;delta;mu];

tau1 = 0.05;
tau2 = 0.05;

noise = [tau1;tau2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial allocation for lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 500;
lambda = zeros(5,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range values for varpi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varpi1 = 0.1;
varpi2 = 0.14;
varpi = linspace(varpi1,varpi2,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Lyapunov exponents lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parfor i = 1:n
    lambda(:,i) = lyapunov_exponents(p,noise,varpi(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

%%%%%%%%%%
h1 = subplot(221);
plot(varpi,lambda(1,:)','r-','LineWidth',2), hold on
plot(varpi,zeros(1,n),'k-','LineWidth',2), hold on
ax = gca;
ax.FontSize = 16;
xlabel('$\varpi$','FontSize',18,'Interpreter','latex'); 
ylabel('$\lambda_1$','FontSize',18,'Interpreter','latex');
ax.TitleHorizontalAlignment = 'left';
title('A','FontSize',20)

%%%%%%%%%%
h2 = subplot(222);
plot(varpi,lambda(2,:)','color',[0 0.5 0],'LineWidth',2), hold on
plot(varpi,zeros(1,n),'k-','LineWidth',2), hold on
ax = gca;
ax.FontSize = 16;
xlabel('$\varpi$','FontSize',18,'Interpreter','latex'); 
ylabel('$\lambda_2$','FontSize',18,'Interpreter','latex');
ax.TitleHorizontalAlignment = 'left';
title('B','FontSize',20)

%%%%%%%%%
h3 = subplot(223);
plot(varpi,lambda(3,:)','color',[0, 0.75, 0.75],'LineWidth',2), hold on
plot(varpi,zeros(1,n),'k-','LineWidth',2), hold on
ax = gca;
ax.FontSize = 16;
xlabel('$\varpi$','FontSize',18,'Interpreter','latex'); 
ylabel('$\lambda_3$','FontSize',18,'Interpreter','latex');
ax.TitleHorizontalAlignment = 'left';
title('C','FontSize',20)

%%%%%%%%%%%%%%%%
h4 = subplot(224);
plot(varpi,lambda(4,:),'b-',varpi,lambda(5,:)','m-','LineWidth',2), hold on
plot(varpi,zeros(1,n),'k-','LineWidth',2), hold on
ax = gca;
ax.FontSize = 16;
xlabel('$\varpi$','FontSize',18,'Interpreter','latex'); 
ylabel('Lyapunov exponent','FontSize',18);
legend('$\lambda_4$','$\lambda_5$','Interpreter','latex','FontSize',13,'Location','best')
ax.TitleHorizontalAlignment = 'left';
title('D','FontSize',20)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','normalized')
set(h4,'position',[.56 .07 .40 .38])
set(h3,'position',[.08 .07 .40 .38])
set(h2,'position',[.56 .55 .40 .38])
set(h1,'position',[.08 .55 .40 .38])
set(gcf, 'PaperSize', [15 10], 'PaperPosition', [0 0 15 10])
print('figure5','-dpdf')