% MAIN FILE TO RUN FIGURE 3 FROM EXAMPLE 2
% -------------------------------------------------------------------
% Clearing memory
close all;
clear all;
clc;
% Local Models of the Fuzzy Subsystems
ri=2;  % Number of local models of the subsystems 
Ns=2;  % Number of subsystems
A = cell(Ns,ri);
A{1,1}=[1 0;0 1];
A{1,2}=[-5 0;0 1];
A{2,2}=[-10 -15;0 -3];
% System parameters
be=1;            % be=beta is a fixed value of \mathcal P
alfa=[0.6 0.4];  % \alpha=[\alpha_1 \alpha_2]
eps1=100;        % \epsilon
eps2=115;
xv=5;            % \x_{\eta}
l=0.49;          % LM= \ell is a real a number defined by designer
mu=.19;          % \mu
fi1=1;
fi2=2;
b=5*sqrt(2);   % b is a real number satisfying |x|\leq b
b1=2;          % b1 is a real number satisfying max_{x\in X}[\sum_{p\in P}\sum_{k\in Gp}h_{pk}]
b2=2;          % b2 is a real number satisfying max_{x\in X}>=(r_p(\sum_{p\in P}\sum_{k\in G_p}h_{pk})+N)
% Calculating the membership function values on the boundary of set Z
h12=(xv^2)/50; h11=1-h12;
h21=0;         h22=1-h21;
% Values used in the grid of Figure 3
alim=[20 35];  % Values for a1
blim=[-10 -5]; % Values for a2

figure(1)
hold on;
for a1=alim(1):(alim(2)-alim(1))/10:alim(2),
    for a2=blim(1):(blim(2)-blim(1))/10:blim(2),
        A{2,1}=[-10 a1;0 a2];
%%     Solving LMIs of Solving LMIs of Valentino et al (2018) (Valentino2018)
        P = lmi_Valentino2018(A,alfa,xv,be,mu,Ns);
        if ~isempty(P)
             plot(a1,a2,'ks','MarkerSize',9)
        end
%%     Solving LMIs of Theorem 2 (eps=100)
         P = Example2_LMI_Theorem2 (A,alfa,fi1,fi2,xv,b,eps1,b1,l,mu,be,b2);
        if ~isempty(P)
              plot(a1,a2,'kx','MarkerSize',9)
        end
%%     Solving LMIs of Theorem 2 (eps=115)
         P = Example2_LMI_Theorem2 (A,alfa,fi1,fi2,xv,b,eps2,b1,l,mu,be,b2);
        if ~isempty(P)
              plot(a1,a2,'k.','MarkerSize',11)
        end
    end
end