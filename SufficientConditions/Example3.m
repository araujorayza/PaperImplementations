% MAIN FILE TO RUN EXAMPLE 3
% -------------------------------------------------------------------
% Clearing memory
clear all;
close all;
clc;
% System parameters
be=1;             % be=beta is a fixed value of \mathcal P
alfa=[0.8 0.2];   % \alpha=[\alpha_1 \alpha_2]
eps=100;          % \epsilon
xv=4;             % \x_{\eta}
mu=.5;            % \mu
xvN=4;            % xvN=2^{r_p}, for p=1 and p=2
L=1.67;           % L is a real number satisfying L<\min_{x\in \partial X}[\sum_{p\in P}\sum_{k\in Gp}h_{pk}]
b=sqrt(32);       % b is a real number satisfying |x|\leq b
b1=1.9;           % b1 is a real number satisfying max_{x\in X}[\sum_{p\in P}\sum_{k\in Gp}h_{pk}]
b2=5.8;           % b2 is a real number satisfying max_{x\in X}>=(r_p(\sum_{p\in P}\sum_{k\in G_p}h_{pk})+N)
LM=L*0.2          % LM= \ell is a real a number defined by designer
fi=[0.2148 0.2148 0.2148 0.2148 0.1616 0.1616 0.1616 0.1616];

% Local Models of the Fuzzy Subsystems
Ns=2; % Number of subsystems
ri=4; % Number of local models of the subsystems 
A=cell(Ns,ri); % Putting the matrices of local models in cell structures

% Matrices of the Local Models
A{1,1}=[5 -1;1 9]; 
A{1,2}=[5 -1;1 -1.87];
A{1,3}=[-1.3 -1;1 9];
A{1,4}=[-1.3 -1;1 -1.87];
A{2,1}=[9 -1;1 19];
A{2,2}=[9 -1;1 -1.77];
A{2,3}=[-8.9 -1;1 16];
A{2,4}=[-8.9 -1;1 -1.87];

celldisp(A) % Displaying the matrices on the screen
% -------------------------------------
% Verifying the feasibility of the set of LMIs
[P,P1,P2] = Example3_LMI (A,alfa,fi,xv,b,LM,eps,b1,b2,mu,be);
if iscell(P)
celldisp(P)
celldisp(P1)
celldisp(P2)   
%% Ploting the level sets of the function V (\Omega_{2\ell} and \Omega_L)
  [x,y] = meshgrid(-4:.05:4,-4:.05:4);
  [tF1,tF2]=size(x);
  for t1=1:tF1,
      for t2=1:tF2,
          w01=(16-x(t1,t2)^2)/(149); w11=1-w01;  
          w02=(32-y(t1,t2)^2)/(198); w12=1-w02;
% Membership Functions (h11=h21=h1; h12=h22=h2;  
%                       h13=h23=h3;  h14=h24=h4)
          h1=w01*w02; h2=w01*w12; h3=w11*w02; h4=w11*w12;
          PH=h2*(P{1,2}+P{2,2})+h4*(P{1,4}+P{2,4});   % PH=h12P12+h14P14+h22P22+h24P24;
          VL(t1,t2) = [x(t1,t2) y(t1,t2)]*PH*[x(t1,t2);y(t1,t2)]; % Function V
      end
  end
  figure(1) % Plotting \Omega_{2\ell} and \Omega_L in Figure 1 
  hold on
  [C1,h]=contour(x,y,VL,[2*LM,2*LM],'r','LineWidth',2);   %\Omega_{2\ell}
  [C2,h]=contour(x,y,VL,[L,L],'k-','LineWidth',2);        % \Omega_L


%% Simulating the nonlinear switched systems with initial condition X01 **

    TF=50;  t=0:.01:TF; % Time of the simulating
    St=length(t);       % Number of elements in the time grid
    X01=[-3;1];         % First initial condition X01
    Taux=[]; Laux=[];
    save Example3.mat Laux Taux;
    save Example3_parameters.mat  A P LM fi
    
    global ip MAX cont % Global variables used in swiched law
    % Starts the simulation from the subsystem (1)
    ip=1;
    
    % Maximum time that a subsystem keeps activated. 
    % The time is count with the iterations
    MAX=30;
    % MAX=200;  % Aparece um ciclo limite
    % MAX=500;  % Aparece um ciclo limite
    
    % Counter to check the maximum time
    cont=0; 
    
    %Finding the numerical solution of the switched system for X01
    y=ode1('Example3_ODE',t,X01);

    % Plotting the phase portrait for X01 
    figure(1)
    hold on
    plot(y(:,1),y(:,2),'b-','linewidth',2)
    plot(X01(1),X01(2),'bo','MarkerSize',8,...
        'MarkerEdgeColor','k','MarkerFaceColor','b')

    % Plotting the switching solution along the time with initial condition X01
    figure(2) 
    hold on
    plot(t,y(:,1),'k--',t,y(:,2),'b-')
    load Example3.mat Laux Taux;
    % Plotting the switching law used to obtain the switching solution    
    plot(Taux,Laux,'b-','linewidth',2)

  % Function V along the switching solution with initial condition X01  
     for i=1:St,
        w01=(16-y(i,1)^2)/(149); w11=1-w01;  
        w02=(32-y(i,2)^2)/(198); w12=1-w02;
        % Membership Functions with constraints
        % h12=h22;    h14=h24;
        h12=w01*w12; h14=w11*w12;
        
        % Function V
        Vc(i)=[y(i,1) y(i,2)]*(h12*(P{1,2}+P{2,2})+...
                        h14*(P{1,4}+P{2,4}))*[y(i,1);y(i,2)];
     end
% Plotting function V along the switching solution with initial condition X01
   figure(3)
   plot(t,Vc,'b','LineWidth',2)
   xlabel('t','fontsize',16)
   ylabel('Vch','fontsize',16)
   axis([-.02 TF 0 0.37])
   
%% Simulating the switched nonlinear systems with initial condition X02 **
     X02=[3;-1.5]; % Second initial condition
     save Example3.mat Laux Taux;
     save Example3_parameters.mat A P LM fi
     ip=1; MAX=50; cont=0; % Setting the global variables
     save Example3.mat Laux Taux;
     St=length(t);
% Finding the numerical solution of the switched system for X02
     y=ode1('Example3_ODE',t,X02);
% Plotting phase portrait for X01 
     figure(1)
     hold on
     plot(y(:,1),y(:,2),'r--','linewidth',2)
     plot(X02(1),X02(2),'ro','MarkerSize',8,...
         'MarkerEdgeColor','k','MarkerFaceColor','r')
% Plotting the switching solution along the time with initial condition X02
     figure(5)
     hold on
     plot(t,y(:,1),'k--',t,y(:,2),'b-')
% Plotting the switching law used to obtain the switching solution
     plot(Taux,Laux,'b','linewidth',2)     
% Function V along the switching solution with initial condition X01  
     for i=1:St,
        w01=(16-y(i,1)^2)/(149); w11=1-w01;  
        w02=(32-y(i,2)^2)/(198); w12=1-w02;
        % Membership Functions with constraints
        % h12=h22;    h14=h24;
        h12=w01*w12; h14=w11*w12;
        
        % Function V
        Vc(i)=[y(i,1) y(i,2)]*(h12*(P{1,2}+P{2,2})+...
                        h14*(P{1,4}+P{2,4}))*[y(i,1);y(i,2)];
     end
% Ploting function V along the switching solution with initial condition X01
   figure(6)
   plot(t,Vc,'b','LineWidth',2)
   xlabel('t','fontsize',16)
   ylabel('Vch','fontsize',16)
   axis([-.02 TF 0 0.37])
 end