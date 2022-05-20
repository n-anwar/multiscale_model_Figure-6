
%% parameters
par.m=.58;      %mosquito per human  (.5838 in White code)
par.a=80/365;     %bitting frequiency (per day)
par.g=1/10;       %mosquito death rate (per day)
par.n=1/12;       %rate of sporogony in mopsquito
par.b=.5;         %transmission probability: mosquito to human
par.c=.23;        %transmission probability: human to mosquito
par.gamma=1/10;       %rate of blood stage infection clearance
par.mu=1/10;  %hypnozoites death rate
par.nu=5;         %Number of hypnozoites per bite
par.alpha=1/1000;  %hypnozoite activation rate
lam=0.03;  %constant lambda (high transmission)
%% discretization
tmax=10*365;      %max simulation time
step=10000;
h=tmax/step;      %step size


imax=100; %max hypnozoite under consideration
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
    
    lambda(j)=lam;
    
    %%Functions within-host model
    
    pA=@(tau) par.alpha*(exp(-par.gamma*tau)-exp(-(par.alpha+par.mu)*tau))/(par.alpha+par.mu-par.gamma);    %eqn (8)
    pH=@(tau) exp(-(par.alpha+par.mu)*tau);  %eqn (7)
    q_t=trapz(t(1:j+1),lambda(1:j+1));          %eqn (11)
    
    NoHyp_int=trapz(t(1:j+1),lambda(1:j+1)./(1+par.nu*exp(-(par.alpha+par.mu)*(t(j+1)-t(1:j+1)))));
    NoHyp=exp(-q_t+NoHyp_int);   %Prob of no hypnozoite ( eqn (14))
    
    NoRel_Inf_int=trapz(t(1:j+1),lambda(1:j+1).*(1-exp(-par.gamma*(t(j+1)-t(1:j+1))))./(1+(par.nu*par.alpha/(par.alpha+par.mu-par.gamma))*(exp(-par.gamma*(t(j+1)-t(1:j+1)))-exp(-(par.alpha+par.mu)*(t(j+1)-t(1:j+1))))));
    NoRel_Inf=exp(-q_t+NoRel_Inf_int);  %Prob of no relapse, no primary infection ( eqn (15))
    
    NoHyp_Rel_Inf_int=trapz(t(1:j+1),lambda(1:j+1).*(1-exp(-par.gamma*(t(j+1)-t(1:j+1))))./(1+(par.nu/(par.alpha+par.mu-par.gamma))*(par.alpha*exp(-par.gamma*(t(j+1)-t(1:j+1)))+(par.mu-par.gamma)*exp(-(par.alpha+par.mu)*(t(j+1)-t(1:j+1))))));
    NoHyp_Rel_Inf=exp(-q_t+NoHyp_Rel_Inf_int);   %Prob of no hypnozoite, no relapse, no primary infection ( eqn (16))
    
    k_T_int=trapz(t(1:j+1),(lambda(1:j+1).*(1-exp(-par.gamma*(t(j+1)-t(1:j+1)))).*(par.nu*pH(t(j+1)-t(1:j+1))))./(1+par.nu*pA(t(j+1)-t(1:j+1))).^2);  %  ( see integration part of eqn (18))
    
    
    Prob_Nohyp(j+1)=NoHyp;
    Prob_NoHyp_Rel_Inf(j+1)=NoHyp_Rel_Inf;
    Prob_NoRel_Inf(j+1)=NoRel_Inf;
    Prob_nohyp_giv_noInf(j+1)=(NoHyp_Rel_Inf)/NoRel_Inf;
    
    
    %% Calculating p (prob of no hypnozoite given infection)
    
    p(j+1)=(NoHyp-NoHyp_Rel_Inf)/(1-NoRel_Inf);  %  (eqn (13))
    
    
    
    Prob_noHyp_giv_No_Inf(j+1)=NoHyp_Rel_Inf/NoRel_Inf;  %Prob of No Hypnozoite given no Infection (P(N_B(t)=0|N_A(t)=N_P(t)=0)
    
    %% Average hypnozoite reservoir size k_T(t)
    
    k_T(j+1)=k_T_int/(1-Prob_noHyp_giv_No_Inf(j+1)); % (see eqn (18))
    %     keyboard
    
    %% Calculating k1(t)
    
    k1_int=trapz(t(1:j+1),(lambda(1:j+1).*(1-exp(-par.gamma*(t(j+1)-t(1:j+1)))).*(par.nu*pH(t(j+1)-t(1:j+1)))./((1+par.nu*(pH(t(j+1)-t(1:j+1))+pA(t(j+1)-t(1:j+1)))).^2))); % eqn (21) for i=1 hypnozoite
    
    g0_g1=NoHyp_Rel_Inf/(NoRel_Inf);  % eqn (21) for i=1 hypnozoite
    
    
    k1(j+1)=g0_g1*k1_int/(1-Prob_noHyp_giv_No_Inf(j+1)); %Prob of 1 hypnozoite given no infection (eqn (20) for i=1)
    
    
    %% definig f(x,y) for ODEs
    dS=@(t,S,I,H) -lambda(j)*S+par.mu*k1(j+1)*H+p(j+1)*par.gamma*I;  %eqn(1)
    dI=@(t,S,I,H) lambda(j)*(S+H)+par.alpha*k_T(j+1)*H-par.gamma*I;  %eqn(2)
    dL=@(t,I,H) -lambda(j)*H-par.mu*k1(j+1)*H-par.alpha*k_T(j+1)*H+(1-p(j+1))*par.gamma*I;  %eqn(3)
    
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
    
    
    l_int_SS(k)=trapz(t(1:dt_step),(lambda(1:dt_step).*(1-exp(-par.gamma*(t(dt_step)-t(1:dt_step)))).*(par.nu*pH(t(dt_step)-t(1:dt_step))).^k)./((1+par.nu*(pH(t(dt_step)-t(1:dt_step))+pA(t(dt_step)-t(1:dt_step)))).^(k+1)));  %(eqn 75 in Mehra et al.  without treatment)
    g_int_SS(k)=trapz(t(1:dt_step),(lambda(1:dt_step).*pH(t(dt_step)-t(1:dt_step)).^k)./(1+par.nu*pH(t(dt_step)-t(1:dt_step))).^(k+1));  % (eqn 79 in Mehra et al. without treatment)
    G_SS(k)=factorial(k)*l_int_SS(k); %reservoir size given no infection  (eqn 78 in Mehra et al.)
    G1_SS(k)=par.nu^k*factorial(k)*g_int_SS(k); %reservoir size  (eqn 74 in Mehra et al.)
    
    
end
BellTerm_SS=IncompleteBellPoly_fac(imax,kmax,G_SS,fac);
BellTerm1_ss=IncompleteBellPoly_fac(imax,kmax,G1_SS,fac);

Last_Bell_SS=zeros(1,imax);
Last_Bell1_ss=zeros(1,imax);
Hyp_L_SS=zeros(1,imax);
Hyp_SS=zeros(1,imax);

for n=1:imax  %loop over max hyp
    for k=1:n
        Last_Bell_SS(k)=BellTerm_SS(n+1,k+1);
        Last_Bell1_ss(k)=BellTerm1_ss(n+1,k+1);
    end
    Hyp_L_SS(n)=sum(Last_Bell_SS,2)*NoHyp_Rel_Inf/(factorial(n)*NoRel_Inf); %Prob of i hyp given no infection
    Hyp_SS(n)=sum(Last_Bell1_ss,2)*NoHyp/factorial(n);  %Prob of hyp i
    
end




% Hyp_I_SS=(Hyp_SS-Hyp_L_SS*NoRel_Inf)/(1-NoRel_Inf);  %prob of i hyp given infection






%% Solution from White model

time=linspace(0,tmax,step+1);
s_h0=zeros(1,imax+1);
i_h0=zeros(1,imax+1);
s_h0(1)=S(1);
i_h0(1)=I(1);
s_m0=S_m(1);
e_m0=E_m(1);
i_m0=I_m(1);

options = odeset('RelTol', 1e-5);

[t_w,y_w]=ode45(@white,[0 tmax],[s_h0,i_h0,s_m0,e_m0,i_m0],options,par,imax,lam);

s_h=y_w(:,1:imax+1);            %Assigning variable
i_h=y_w(:,imax+2:2*imax+2);

I_white=sum(i_h,2);
L_white=sum(s_h(:,2:end),2);


%% Figure

subplot(1,2,1)
plot(t,I,'b',t,L,'--b','linewidth',2)
hold on
plot(t_w,I_white,'r',t_w,L_white','--r','linewidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Time (Days)','fontweight','bold','Fontsize',20)
ylabel('Fraction of human population','fontweight','bold','Fontsize',20)
legend('I from our multiscale model','L from our multiscale model','I from 2(L_{max}+1) model','L from from 2(L_{max}+1) model','fontsize',17)
xlim([0 100])
ylim([0 .059])
title('A')

%% prealocating
Hyp_dist_SIL=zeros(imax+1,1); %probability of having hypnozoite from our model
Hyp_dist_white=zeros(imax+1,1); %probability of having hypnozoite from our model

Hyp_dist=zeros(imax+1,2); %prealocation for bar plot (side by side)

%% finding probability
Hyp_dist_SIL(1)=Prob_Nohyp(end); %probability of having 0 hypnozoite at steady state

for i=1:imax
    Hyp_dist_SIL(i+1)=Hyp_SS(i);  %probability of having i hypnozoite at steady state
end

Hyp_dist_white(1)=s_h(end,1)+i_h(end,1);    %probability of having 0 hypnozoite from white model


Hyp_dist(1,:)=[Hyp_dist_SIL(1) Hyp_dist_white(1)]; %probability of having 0 hypnozoite from both model

for i=1:imax
    Hyp_dist_white(i+1)=s_h(end,i+1)+i_h(end,i+1); %probability of having i hypnozoite from white model
    Hyp_dist(i+1,:)=[Hyp_dist_SIL(i+1) Hyp_dist_white(i+1)];  %probability of having i hypnozoite from both model
end




subplot(1,2,2)
hyp_number=0:imax;
bar(hyp_number,Hyp_dist);
ax = gca;
ax.FontSize = 15;
xlabel('Number of Hypnozoites within individuals','fontweight','bold','Fontsize',20)
ylabel('Probability','fontweight','bold','Fontsize',20)
legend('Our multi-scale model','2(L_{max}+1) model','fontsize',17)
xlim([-1 21])
title('B')



%% ODE function of White model
function ydot=white(~,y,par,imax,lam)

ds_h=zeros(imax+1,1);
di_h=zeros(imax+1,1);


%geometric distribution of number of hypnozoites per infection
mat=zeros(imax+1,imax+1);
prob=1/(par.nu+1);
for i=1:imax
    for j=1:i
        mat(i,j)=prob*(1-prob).^(i-j);
        
    end
end

mat(imax+1,:)=1-sum(mat,1); % for balanceing the total population (1-column sum(all previous))


s_h=y(1:imax+1);            %Assigning variable
i_h=y(imax+2:2*imax+2);

s_m=y(2*imax+3);
e_m=y(2*imax+4);
i_m=y(2*imax+5);

lambda=lam;
lambda_mat=mat*lambda;

if imax==0
    s_h(2)=0;
    i_h(2)=0;
end

ds_h(1)=-lambda*s_h(1)+par.gamma*i_h(1)+par.mu*s_h(2);
di_h(1)=-lambda*i_h(1)-par.gamma*i_h(1)+(par.mu+par.alpha)*i_h(2)+par.alpha*s_h(2)...
    +(lambda_mat(1,1)*(s_h(1)+i_h(1)));

for i=2:imax
    
    ds_h(i)=-lambda*s_h(i)+par.gamma*i_h(i)-(i-1)*(par.mu+par.alpha)*s_h(i)+i*par.mu*s_h(i+1);
    di_h(i)=-lambda*i_h(i)-par.gamma*i_h(i)-(i-1)*(par.mu+par.alpha)*i_h(i)+i*(par.mu+par.alpha)*i_h(i+1)+i*par.alpha*s_h(i+1)...
        +sum(lambda_mat(i,:)*(s_h+i_h));
    
end

ds_h(imax+1)=-lambda*s_h(imax+1)+par.gamma*i_h(imax+1)-imax*(par.mu+par.alpha)*s_h(imax+1);
di_h(imax+1)=-lambda*i_h(imax+1)-par.gamma*i_h(imax+1)-imax*(par.mu+par.alpha)*i_h(imax+1)...
    +sum(lambda_mat(imax+1,:)*(s_h+i_h));

ds_m=par.g-par.a*par.c*sum(i_h)*s_m-par.g*s_m;
de_m=par.a*par.c*sum(i_h)*s_m-e_m*(1/par.n)-par.g*e_m;
di_m=e_m*(1/par.n)-par.g*i_m;


ydot=[ds_h;di_h;ds_m;de_m;di_m];


end