
%% parameters
par.m=.5838; %mosquito per human  (.5838 in White code)
par.a=80/365;% %; %bitting frequiency (per day)
par.g=1/10;      %mosquito death rate (per day)
par.n=1/12;        %rate of sporogony in mopsquito
par.b=.5;        %transmission probability: mosquito to human
par.c=.23;       %transmission probability: human to mosquito
par.r=1/60;        %rate of blood stage infection clearance
par.omega=1/425  ;   %hypnozoites death rate
par.nu=5;   %Number of hypnozoites per infection
par.alpha=1/332;%1/332;

%% discretization
tmax=10*365;
step=150000;
h=tmax/step;   %step size


imax=100;       %Bell parameter
kmax=imax+5;    %Bell parameter
%% Initials
S=zeros(step+1,1);
I=zeros(step+1,1);
L=zeros(step+1,1);
S_m=zeros(step+1,1);
I_m=zeros(step+1,1);
E_m=zeros(step+1,1);

S_m(1)=.95;
E_m(1)=0;
I_m(1)=.05;

t=zeros(step+1,1);

S(1)=1;
I(1)=0;
L(1)=0;



%% Pre-alocating space
G_SS=zeros(1,kmax);
G1_SS=zeros(1,kmax);
% G=zeros(1,kmax);
g_int_SS=zeros(1,kmax);
g1_int_SS=zeros(1,kmax);
Bell=zeros(step+1,1);

Prob_Nohyp=zeros(step+1,1);
Prob_NoHyp_Rel_Inf=zeros(step+1,1);
Prob_NoRel_Inf=zeros(step+1,1);
p=zeros(step+1,1);

lambda=zeros(step+1,1);

s0=zeros(step+1,1);
Tot_prob=zeros(step+1,1);

Bell_SS=zeros(1,imax);
Bell1_SS=zeros(1,imax);
k1=zeros(step+1,1);
T_sum=zeros(step+1,1);

Prob_nohyp_giv_noInf=zeros(step+1,1);
fac=zeros(1,imax+1);
fac(1)=1;                   %pre calculating factorial for Bell_poly function that will reduce computation time
for j=1:imax
    fac(j+1)=factorial(j);
end


%% RK4 algorithm
for j=1:step
    t(j+1)=j*h;
    
    lambda(j)=par.m*par.a*par.b*I_m(j);
    
    %%Functions within-host model
    pA=@(tau) par.alpha*(exp(-par.r*tau)-exp(-(par.alpha+par.omega)*tau))/(par.alpha+par.omega-par.r);    %eqn (12)
    pH=@(tau) exp(-(par.alpha+par.omega)*tau);  %eqn (11)
    q_t=trapz(t(1:j+1),lambda(1:j+1));
    
    NoHyp_int=trapz(t(1:j+1),lambda(1:j+1)./(1+par.nu*exp(-(par.alpha+par.omega)*(t(j+1)-t(1:j+1))))); 
    NoHyp=exp(-q_t+NoHyp_int);   %Prob of no hypnozoite ( see eqn (17))
 
    NoRel_Inf_int=trapz(t(1:j+1),lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1))))./(1+(par.nu*par.alpha/(par.alpha+par.omega-par.r))*(exp(-par.r*(t(j+1)-t(1:j+1)))-exp(-(par.alpha+par.omega)*(t(j+1)-t(1:j+1))))));
    NoRel_Inf=exp(-q_t+NoRel_Inf_int);
    
    NoHyp_Rel_Inf_int=trapz(t(1:j+1),lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1))))./(1+(par.nu/(par.alpha+par.omega-par.r))*(par.alpha*exp(-par.r*(t(j+1)-t(1:j+1)))+(par.omega-par.r)*exp(-(par.alpha+par.omega)*(t(j+1)-t(1:j+1))))));
    NoHyp_Rel_Inf=exp(-q_t+NoHyp_Rel_Inf_int);
    
    T_sum_int=trapz(t(1:j+1),(lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1)))).*(par.nu*pH(t(j+1)-t(1:j+1))))./(1+par.nu*pA(t(j+1)-t(1:j+1))).^2);%  ( see eqn (24))
    
    
    Prob_Nohyp(j+1)=NoHyp;
    %Prob of No Relapse and no Infection ( see eqn (18))
    
    %Prob of No Hypnozoite, no Relapse and no Infection ( see eqn (19))
    
    Prob_NoHyp_Rel_Inf(j+1)=NoHyp_Rel_Inf;
    Prob_NoRel_Inf(j+1)=NoRel_Inf;
    Prob_nohyp_giv_noInf(j+1)=(NoHyp_Rel_Inf)/NoRel_Inf;
    
    
    %% Calculating p (prob of no hypnozoite given infection
    if NoRel_Inf==1
        p(j+1)=1;
    else
        p(j+1)=(NoHyp-NoHyp_Rel_Inf)/(1-NoRel_Inf);  %  (eqn (8))
    end
    
    
    s0(j+1)=NoHyp_Rel_Inf/NoRel_Inf;  %Prob of No Hypnozoite given no Infection (P(N_B(t)=0|N_A(t)=N_P(t)=0) (see eqn 20)
    
    %% (k1+2k2+3k3+...upto infty)
    if s0(j+1)==1
        T_sum(j+1)=0;
    else
        T_sum(j+1)=T_sum_int/(1-s0(j+1)); %k1+2k2+3k3+...+ipto infty (see eqn (24))
    end
    
    Tot_prob(j+1)=T_sum(j+1);
    
    
    % Calculating k1(t)
    
        g_int=trapz(t(1:j+1),(lambda(1:j+1).*(1-exp(-par.r*(t(j+1)-t(1:j+1)))).*(par.nu*pH(t(j+1)-t(1:j+1)))./((1+par.nu*(pH(t(j+1)-t(1:j+1))+pA(t(j+1)-t(1:j+1)))).^2)));
    G=g_int;   %eqn 15
    
    
    BellTerm=IncompleteBellPoly_fac(1,1,G,fac);   % calculating bell polynomial where we only need the last row.
    
    
    Last_Bell=BellTerm(2,2); %only need the last row
    
    Bell(j+1)=sum(Last_Bell,2)*NoHyp_Rel_Inf/(NoRel_Inf);  % eqn (21) for i=1 hypnozoite
    
    if s0(j+1)==1
        k1(j+1)=0;
    else
        k1(j+1)=Bell(j+1)/(1-s0(j+1)); %Prob of 1 hypnozoite given no infection (eqn (20) for i=1)
    end
    
    %% definig f(x,y) for ODEs
    dS=@(t,S,I,H) -lambda(j)*S+par.omega*k1(j+1)*H+p(j+1)*par.r*I;  %eqn(25)
    dI=@(t,S,I,H) lambda(j)*(S+H)+par.alpha*Tot_prob(j+1)*H-par.r*I;  %eqn(26)
    dL=@(t,I,H) -lambda(j)*H-par.omega*k1(j+1)*H-par.alpha*Tot_prob(j+1)*H+(1-p(j+1))*par.r*I;  %eqn(27)
    
    dS_m=@(t,S_m,I) par.g-par.a*par.c*I*S_m-par.g*S_m;  %eqn(28)
    dE_m=@(t,S_m,I,E_m) par.a*par.c*I*S_m-(par.g+par.n)*E_m;  %eqn(29)
    dI_m=@(t,I_m,E_m) par.n*E_m-par.g*I_m; %eqn(30)
    
    
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