
%% parameters
par.a=80/365;     %bitting frequiency (per day)
par.g=1/10;       %mosquito death rate (per day)
par.r=1/60;       %rate of blood stage infection clearance
par.omega=1/425;  %hypnozoites death rate
par.nu=5;         %Number of hypnozoites per bite
par.alpha=1/332;  %hypnozoites death rate

%% discretization 
tmax=1000;  %total simulation time
step=5000;   %number of steps
h=tmax/step;   %step size

%% Pre-alocating space

t=zeros(step,1);
p=zeros(step+1,1);                    %probability of no hypnozoites given blood-stage infection

lambda=zeros(step+1,1);               %Force of reinfection

Prob_nohyp_giv_noInf=zeros(step+1,1);

k_T=zeros(step+1,1);                  %Averaze size of hypnozoite reservoir

k1=zeros(step+1,1);                   %probability of having 1 hypnozoite given no blood-stage infection


for j=1:step
t(j+1)=j*h;
    
    lambda(j)=0.005;
    
    %%Functions within-host model
    
    pH=@(tau) exp(-(par.alpha+par.omega)*tau);  %eqn (7)
    pA=@(tau) par.alpha*(exp(-par.r*tau)-exp(-(par.alpha+par.omega)*tau))/(par.alpha+par.omega-par.r);    %eqn (8)
    
    pC=@(tau) (par.alpha/(par.alpha+par.omega))*(1-pH(tau))-pA(tau); %eqn (9)
    
    pD=@(tau) (par.omega/(par.alpha+par.omega))*(1-pH(tau));%eqn (10)
  
    q_t=trapz(t(1:j+1),lambda(1:j+1));          %eqn (11)
    
    NoHyp_int=trapz(t(1:j+1),lambda(1:j+1)./(1+par.nu*exp(-(par.alpha+par.omega)*(t(j+1)-t(1:j+1)))));
    NoHyp=exp(-q_t+NoHyp_int);   %Prob of no hypnozoite ( eqn (14))
    
    NoRel_Inf_int=trapz(t(1:j+1),lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1))))./(1+(par.nu*par.alpha/(par.alpha+par.omega-par.r))*(exp(-par.r*(t(j+1)-t(1:j+1)))-exp(-(par.alpha+par.omega)*(t(j+1)-t(1:j+1))))));
    NoRel_Inf=exp(-q_t+NoRel_Inf_int);  %Prob of no relapse, no primary infection ( eqn (15))
    
    NoHyp_Rel_Inf_int=trapz(t(1:j+1),lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1))))./(1+(par.nu/(par.alpha+par.omega-par.r))*(par.alpha*exp(-par.r*(t(j+1)-t(1:j+1)))+(par.omega-par.r)*exp(-(par.alpha+par.omega)*(t(j+1)-t(1:j+1))))));
    NoHyp_Rel_Inf=exp(-q_t+NoHyp_Rel_Inf_int);   %Prob of no hypnozoite, no relapse, no primary infection ( eqn (16))
    
    k_T_int=trapz(t(1:j+1),(lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1)))).*(par.nu*pH(t(j+1)-t(1:j+1))))./(1+par.nu*pA(t(j+1)-t(1:j+1))).^2);  %  ( see integration part of eqn (18))
 
    Prob_nohyp_giv_noInf(j+1)=(NoHyp_Rel_Inf)/NoRel_Inf;
    
 %% Calculating p (prob of no hypnozoite given infection)
 
    p(j+1)=(NoHyp-NoHyp_Rel_Inf)/(1-NoRel_Inf);  %  (eqn (13))
    
     %% Average hypnozoite reservoir size k_T(t)
    
    k_T(j+1)=k_T_int/(1-Prob_nohyp_giv_noInf(j+1)); % (see eqn (18))

     %% Calculating k1(t)
    
    k1_int=trapz(t(1:j+1),(lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1)))).*(par.nu*pH(t(j+1)-t(1:j+1)))./((1+par.nu*(pH(t(j+1)-t(1:j+1))+pA(t(j+1)-t(1:j+1)))).^2))); % integral part of eqn (17) 
    
    g0_g1=NoHyp_Rel_Inf/(NoRel_Inf);  % eqn (17)
    
    
    k1(j+1)=g0_g1*k1_int/(1-Prob_nohyp_giv_noInf(j+1)); %Prob of 1 hypnozoite given no infection (eqn (20) for i=1)
    
    
end

%% Figure 4
% figure
% subplot(1,2,1)
plot(t,pH(t),t,pA(t),t,pD(t),t,pC(t),'linewidth',2)
title('A')
ax = gca;
ax.FontSize = 18; 
xlabel('Time (Days)','fontweight','bold','Fontsize',22)
ylabel('Probability','fontweight','bold','Fontsize',22)
hl = legend('$p_B(t)$','$p_A(t)$','$p_D(t)$','$p_C(t)$','fontsize',22);
set(hl, 'Interpreter','latex')
hold on
% subplot(1,2,2)
% yyaxis left
% plot(t(2:end),p(2:end),'linewidth',2)
% hold on
% plot(t(2:end),k1(2:end),'linewidth',2)
% 
% ylabel('Probability','fontweight','bold','Fontsize',22)
% % hold on
% yyaxis right
% plot(t(2:end),k_T(2:end),'linewidth',2)
% title('B')
% 
% ax = gca;
% ax.FontSize = 18;
% xlabel('Time (days)','fontweight','bold','Fontsize',22)
% ylabel('$k_T$','fontweight','bold','Fontsize',22,'Interpreter','latex')
% % ylabel('$Total relapse rate$','fontweight','bold','Fontsize',20,'Interpreter','Latex')
% hl = legend('$p(t)$','$k_1(t)$','$ k_T(t)$','fontsize',22);
% set(hl, 'Interpreter','latex')

