
%% parameters
par.m=.5838;      %mosquito per human  (.5838 in White code)
par.a=80/365;     %bitting frequiency (per day)
par.g=1/10;       %mosquito death rate (per day)
par.n=1/12;       %rate of sporogony in mopsquito
par.b=.5;         %transmission probability: mosquito to human
par.c=.23;        %transmission probability: human to mosquito
par.r=1/60;       %rate of blood stage infection clearance
par.omega=1/425;  %hypnozoites death rate
par.nu=5;         %Number of hypnozoites per infection
par.alpha=1/332;  %hypnozoites death rate

%% discretization
tmax=10*365;      %max simulation time
step=100000;
h=tmax/step;      %step size

%% Pre-alocating space

t=zeros(step+1,1);
S=zeros(step+1,1);                    %fraction of susceptibles
I=zeros(step+1,1);                    %fraction of blood-stage infected
L=zeros(step+1,1);                    %fraction of liver-stage infected
S_m=zeros(step+1,1);                  %fraction of susceptibles mosquitoes
I_m=zeros(step+1,1);                  %fraction of infected mosquitoes
E_m=zeros(step+1,1);                  %fraction of exposed mosquitoes

Prob_Nohyp=zeros(step+1,1);           %probability of no hypnozoites
Prob_NoHyp_Rel_Inf=zeros(step+1,1);   %probability of no hypnozoites, no relapse, and no primary infection
Prob_NoRel_Inf=zeros(step+1,1);       %probability of no relapse and no  primary infection
p=zeros(step+1,1);                    %probability of no hypnozoites given blood-stage infection

lambda=zeros(step+1,1);               %Force of reinfection

Prob_noHyp_giv_No_Inf=zeros(step+1,1);
k_T=zeros(step+1,1);                  %Averaze size of hypnozoite reservoir

k1=zeros(step+1,1);                   %probability of having 1 hypnozoite given no blood-stage infection

Prob_nohyp_giv_noInf=zeros(step+1,1); %probability of having no hypnozoite given no blood-stage infection

%% Initials conditions

S_m(1)=.95;
E_m(1)=0;
I_m(1)=.05;
S(1)=1;
I(1)=0;
L(1)=0;
%% RK4 algorithm
for j=1:step
    t(j+1)=j*h;
    
    lambda(j)=par.m*par.a*par.b*I_m(j);
    
    %%Functions within-host model
    
    pA=@(tau) par.alpha*(exp(-par.r*tau)-exp(-(par.alpha+par.omega)*tau))/(par.alpha+par.omega-par.r);    %eqn (8)
    pH=@(tau) exp(-(par.alpha+par.omega)*tau);  %eqn (7)
    q_t=trapz(t(1:j+1),lambda(1:j+1));          %eqn (11)
    
    NoHyp_int=trapz(t(1:j+1),lambda(1:j+1)./(1+par.nu*exp(-(par.alpha+par.omega)*(t(j+1)-t(1:j+1)))));
    NoHyp=exp(-q_t+NoHyp_int);   %Prob of no hypnozoite ( eqn (14))
    
    NoRel_Inf_int=trapz(t(1:j+1),lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1))))./(1+(par.nu*par.alpha/(par.alpha+par.omega-par.r))*(exp(-par.r*(t(j+1)-t(1:j+1)))-exp(-(par.alpha+par.omega)*(t(j+1)-t(1:j+1))))));
    NoRel_Inf=exp(-q_t+NoRel_Inf_int);  %Prob of no relapse, no primary infection ( eqn (15))
    
    NoHyp_Rel_Inf_int=trapz(t(1:j+1),lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1))))./(1+(par.nu/(par.alpha+par.omega-par.r))*(par.alpha*exp(-par.r*(t(j+1)-t(1:j+1)))+(par.omega-par.r)*exp(-(par.alpha+par.omega)*(t(j+1)-t(1:j+1))))));
    NoHyp_Rel_Inf=exp(-q_t+NoHyp_Rel_Inf_int);   %Prob of no hypnozoite, no relapse, no primary infection ( eqn (16))
    
    k_T_int=trapz(t(1:j+1),(lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1)))).*(par.nu*pH(t(j+1)-t(1:j+1))))./(1+par.nu*pA(t(j+1)-t(1:j+1))).^2);  %  ( see integration part of eqn (18))
    
    
    Prob_Nohyp(j+1)=NoHyp;
    Prob_NoHyp_Rel_Inf(j+1)=NoHyp_Rel_Inf;
    Prob_NoRel_Inf(j+1)=NoRel_Inf;
    Prob_nohyp_giv_noInf(j+1)=(NoHyp_Rel_Inf)/NoRel_Inf;
    
    
    %% Calculating p (prob of no hypnozoite given infection)
    
    p(j+1)=(NoHyp-NoHyp_Rel_Inf)/(1-NoRel_Inf);  %  (eqn (13))
    
    
    
    Prob_noHyp_giv_No_Inf(j+1)=NoHyp_Rel_Inf/NoRel_Inf;  %Prob of No Hypnozoite given no Infection (P(N_B(t)=0|N_A(t)=N_P(t)=0)
    
    %% Average hypnozoite reservoir size k_T(t)
    
    k_T(j+1)=k_T_int/(1-Prob_noHyp_giv_No_Inf(j+1)); % (see eqn (18))
    
    
    %% Calculating k1(t)
    
    k1_int=trapz(t(1:j+1),(lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1)))).*(par.nu*pH(t(j+1)-t(1:j+1)))./((1+par.nu*(pH(t(j+1)-t(1:j+1))+pA(t(j+1)-t(1:j+1)))).^2)));
    
    g0_g1=NoHyp_Rel_Inf/(NoRel_Inf);  % eqn (21) for i=1 hypnozoite
    
    
    k1(j+1)=g0_g1*k1_int/(1-Prob_noHyp_giv_No_Inf(j+1)); %Prob of 1 hypnozoite given no infection (eqn (20) for i=1)
    
    
    %% definig f(x,y) for ODEs
    dS=@(t,S,I,H) -lambda(j)*S+par.omega*k1(j+1)*H+p(j+1)*par.r*I;  %eqn(1)
    dI=@(t,S,I,H) lambda(j)*(S+H)+par.alpha*k_T(j+1)*H-par.r*I;  %eqn(2)
    dL=@(t,I,H) -lambda(j)*H-par.omega*k1(j+1)*H-par.alpha*k_T(j+1)*H+(1-p(j+1))*par.r*I;  %eqn(3)
    
    dS_m=@(t,S_m,I) par.g-par.a*par.c*I*S_m-par.g*S_m;  %eqn(4)
    dE_m=@(t,S_m,I,E_m) par.a*par.c*I*S_m-(par.g+par.n)*E_m;  %eqn(5)
    dI_m=@(t,I_m,E_m) par.n*E_m-par.g*I_m; %eqn(6)
    
    
    %% RK4 method
    
    Sk1=dS(t(j),S(j),I(j),L(j));
    Ik1=dI(t(j),S(j),I(j),L(j));
    Lk1=dL(t(j),I(j),L(j));
    Smk1=dS_m(t(j),S_m(j),I(j));
    Emk1=dE_m(t(j),S_m(j),I(j),E_m(j));
    Imk1=dI_m(t(j),I_m(j),E_m(j));
    
    Sk2=dS(t(j)+.5*h,S(j)+Sk1*.5*h,I(j)+Ik1*.5*h,L(j)+Lk1*.5*h);
    Ik2=dI(t(j)+.5*h,S(j)+Sk1*.5*h,I(j)+Ik1*.5*h,L(j)+Lk1*.5*h);
    Lk2=dL(t(j)+.5*h,I(j)+Ik1*.5*h,L(j)+Lk1*.5*h);
    Smk2=dS_m(t(j)+.5*h,S_m(j)+Smk1*.5*h,I(j)+Ik1*.5*h);
    Emk2=dE_m(t(j)+.5*h,S_m(j)+Smk1*.5*h,I(j)+Ik1*.5*h,E_m(j)+Emk1*.5*h);
    Imk2=dI_m(t(j)+.5*h,I_m(j)+Imk1*.5*h,E_m(j)+Emk1*.5*h);
    
    
    Sk3=dS(t(j)+.5*h,S(j)+Sk2*.5*h,I(j)+Ik2*.5*h,L(j)+Lk2*.5*h);
    Ik3=dI(t(j)+.5*h,S(j)+Sk2*.5*h,I(j)+Ik2*.5*h,L(j)+Lk2*.5*h);
    Lk3=dL(t(j)+.5*h,I(j)+Ik2*.5*h,L(j)+Lk2*.5*h);
    Smk3=dS_m(t(j)+.5*h,S_m(j)+Smk2*.5*h,I(j)+Ik2*.5*h);
    Emk3=dE_m(t(j)+.5*h,S_m(j)+Smk2*.5*h,I(j)+Ik2*.5*h,E_m(j)+Emk2*.5*h);
    Imk3=dI_m(t(j)+.5*h,I_m(j)+Imk2*.5*h,E_m(j)+Emk2*.5*h);
    
    
    Sk4=dS(t(j)+h,S(j)+Sk3*h,I(j)+Ik3*h,L(j)+Lk3*h);
    Ik4=dI(t(j)+h,S(j)+Sk3*h,I(j)+Ik3*h,L(j)+Lk3*h);
    Lk4=dL(t(j)+h,I(j)+Ik3*h,L(j)+Lk3*h);
    Smk4=dS_m(t(j)+h,S_m(j)+Smk3*h,I(j)+Ik3*h);
    Emk4=dE_m(t(j)+h,S_m(j)+Smk3*h,I(j)+Ik3*h,E_m(j)+Emk3*h);
    Imk4=dI_m(t(j)+h,I_m(j)+Imk1*.5*h,E_m(j)+Emk3*h);
    
    
    S(j+1)=S(j)+h*(Sk1+2*Sk2+2*Sk3+Sk4)/6;
    I(j+1)=I(j)+h*(Ik1+2*Ik2+2*Ik3+Ik4)/6;
    L(j+1)=L(j)+h*(Lk1+2*Lk2+2*Lk3+Lk4)/6;
    S_m(j+1)=S_m(j)+h*(Smk1+2*Smk2+2*Smk3+Smk4)/6;
    E_m(j+1)=E_m(j)+h*(Emk1+2*Emk2+2*Emk3+Emk4)/6;
    I_m(j+1)=I_m(j)+h*(Imk1+2*Imk2+2*Imk3+Imk4)/6;
    
    
end

%% Hypnozoite distribution
dist_time=t(end);
dt_step=round(dist_time/h);
imax=100;       %Bell parameter
kmax=imax+5;    %Bell parameter
G_SS=zeros(1,kmax);
G1_SS=zeros(1,kmax);
g_int_SS=zeros(1,kmax);
l_int_SS=zeros(1,kmax);

fac=zeros(1,imax+1);
fac(1)=1;                   %pre calculating factorial for Bell_poly function that will reduce computation time
for j=1:imax
    fac(j+1)=factorial(j);
end

for k=1:kmax
    %% hyp dist in SS  g for P(H=i|A=P=0), g1 for P(H=i), g2=P(H=i|P=0)
    
    
    l_int_SS(k)=trapz(t(1:dt_step),(lambda(1:dt_step).*(1-exp(-par.r*(t(dt_step)-t(1:dt_step)))).*(par.nu*pB(t(dt_step)-t(1:dt_step))).^k)./((1+par.nu*(pB(t(dt_step)-t(1:dt_step))+pA(t(dt_step)-t(1:dt_step)))).^(k+1)));  %(eqn 75 in Mehra et al.  without treatment)
    g_int_SS(k)=trapz(t(1:dt_step),(lambda(1:dt_step).*pB(t(dt_step)-t(1:dt_step)).^k)./(1+par.nu*pB(t(dt_step)-t(1:dt_step))).^(k+1));  % (eqn 79 in Mehra et al. without treatment)
    G_SS(k)=factorial(k)*l_int_SS(k); %reservoir size given no infection  (eqn 78 in Mehra et al.)
    G1_SS(k)=par.nu^k*factorial(k)*g_int_SS(k); %reservoir size  (eqn 74 in Mehra et al.)
    
    
end
BellTerm_SS=IncompleteBellPoly_fac(imax,kmax,G_SS,fac);
BellTerm1_ss=IncompleteBellPoly_fac(imax,kmax,G1_SS,fac);

Last_Bell_SS=zeros(1,imax);
Last_Bell1_ss=zeros(1,imax);
Bell_SS=zeros(1,imax);
Bell1_SS=zeros(1,imax);

for n=1:imax  %loop over max hyp
    for k=1:n
        Last_Bell_SS(k)=BellTerm_SS(n+1,k+1);
        Last_Bell1_ss(k)=BellTerm1_ss(n+1,k+1);
        %             Last_Bell2_ss(l)=BellTerm2_ss(n+1,l+1);
    end
    Bell_SS(n)=sum(Last_Bell_SS,2)*NoHyp_Rel_Inf/(factorial(n)*NoRel_Inf); %Prob of i hyp given no infection
    Bell1_SS(n)=sum(Last_Bell1_ss,2)*NoHyp/factorial(n);  %Prob of hyp i
    %         Bell2_noBI(i+1,n)=sum(Last_Bell2_ss,2)*NoHyp_Inf/(NoInf*factorial(n));   %Prob of hyp i |no primary infection
    
end




Hyp_dist_giv_infection=(Bell1_SS-Bell_SS*NoRel_Inf)/(1-NoRel_Inf);  %prob of i hyp given infection



Hyp_dist=zeros(imax+1,1);

Hyp_dist(1)=NoHyp;

for i=1:imax
    Hyp_dist(i+1)=Bell1_SS(i);
end


%% Figure 6
subplot(2,2,1) %Figure 6A

hold on
plot(t,I,'linewidth',2)
plot(t,L,'linewidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Time (Days)','fontweight','bold','Fontsize',20)
ylabel('Fractions','fontweight','bold','Fontsize',20)
legend('Blood-stsge infected (I)','Liver-stsge infected (L)')
xlim([0 tmax])
title('A')



subplot(2,2,2)  %Figure 6B

Hyp_dist=zeros(imax+1,1);
Hyp_dist(1)=Prob_Nohyp(end);
for i=1:imax
    Hyp_dist(i+1)=Bell1_SS(i);
end
hyp_number=0:imax;
bar(hyp_number(1:end),Hyp_dist(1:end));
ax = gca;
ax.FontSize = 15;
xlabel('Number of Hypnozoites','fontweight','bold','Fontsize',20)
ylabel('Probability','fontweight','bold','Fontsize',20)
xlim([-1 40])
title('B')



subplot(2,2,3)   %Figure 6C

Hyp_dist_H=Bell_SS(:)/(1-Prob_nohyp_giv_noInf(end));
hyp_number_H=1:imax;
bar(hyp_number_H(1:end),Hyp_dist_H(1:end));
ax = gca;
ax.FontSize = 15;
xlabel('Hypnozoites in liver-stage infected individuals','fontweight','bold','Fontsize',20)
ylabel('Probability','fontweight','bold','Fontsize',20)
xlim([0 40])
title('C')



subplot(2,2,4) %Figure 6D

Hyp_dist_I=zeros(imax+1,1);
Hyp_dist_I(1)=p(end);
Hyp_dist_I(2:end)=Hyp_dist_giv_infection(:);
bar(hyp_number(1:end),Hyp_dist_I(1:end));
ax = gca;
ax.FontSize = 15;
xlabel('Hypnozoites in blood-stage infected individuals','fontweight','bold','Fontsize',20)
ylabel('Probability','fontweight','bold','Fontsize',20)
xlim([-1 40])
title('D')
