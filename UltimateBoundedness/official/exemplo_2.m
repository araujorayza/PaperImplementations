
% ROTINA USADA PARA SIMULAR O SISTEMA FUZZY DO EXEMPLO 1
% -------------------------------------------------------------------
% Limpeza de memória
close all;
clear all;
clc;
% PARÂMETROS DO PROBLEMA
ri=2;
Ns=2;
% Parâmetros do sistema
% Parametros usados na combinação convexa
alfa=[0.6 0.4];
beta=1;
xv=5;
gama=0.19;  % Usado nas LMIs para limitar as matrizes P_pk
lamb=.02; % Restriçao na taxa de decaimento
A = cell(Ns,ri);

% Modelo ALTERADO 10/03/2017
A{1,1}=[1 0;0 1];
A{1,2}=[-5 0;0 1];
A{2,1}=[-10 15;0 -3];
A{2,2}=[-10 -15;0 -3];

% % Modelo ALTERADO 08/04/2017
% A{1,1}=[1 0;0 1];
% A{1,2}=[-5 -10;0 1];
% A{2,1}=[-10 15;0 -3];
% A{2,2}=[-10 -15;0 -3];

% -------------------------------------------------------------------
P = lmi_chav_fuzzychaveado(A,alfa,beta,xv,gama);
% P = lmi_chav_fuzzychaveado2gama (A,alfa,beta,xv,gama);
if ~isempty(P)
  celldisp(A)
  celldisp(P)
%   AutoP11=eig(P{1,1})
%   AutoP12=eig(P{1,2})
%   AutoP21=eig(P{2,1})
%   AutoP22=eig(P{2,2})
  dlim = [-xv xv];
  % Exibindo a região onde Vponto é positiva e calculando o máximo da
  % da V nessa região
  [x1,x2] = meshgrid(dlim(1):.1:dlim(2),dlim(1):.1:dlim(2));
  [tF1,tF2]=size(x1);
  figure(1);
  hold on;
  VO=[];  DVH=[];
  for t1=1:tF1,
     for t2=1:tF2,
        XT=[x1(t1,t2);x2(t1,t2)]; % Vetor de estado.
        % Funções de pertinência
        h12=(x1(t1,t2)^2+x2(t1,t2)^2)/50; h11=1-h12;
        h21=(x2(t1,t2)^2)/25;  h22=1-h21;

        PT = h12*P{1,2}+h21*P{2,1};
        VO(t1,t2) =  XT'*PT*XT;
        Ah1=h11*A{1,1}+h12*A{1,2};
        Ah2=h21*A{2,1}+h22*A{2,2};        
        S1=Ah1*XT; % Subsistema 1
        S2=Ah2*XT; % Subsistema 2
         S=alfa(1)*S1+alfa(2)*S2;
        dh12=(x1(t1,t2)*(alfa(1)*S1(1)+alfa(2)*S2(1))+x2(t1,t2)*(alfa(1)*S1(2)+alfa(2)*S2(2)))/25; % h12ponto
        dh21=2*x2(t1,t2)*(alfa(1)*S1(2)+alfa(2)*S2(2))/25;
        PHX=XT'*(dh12*P{1,2}+dh21*P{2,1})*XT;
        if  (PHX > 0)
           %plot(x1(t1,t2),x2(t1,t2),'b*','MarkerSize',3)
           plot(x1(t1,t2),x2(t1,t2),'*','Color',[0.6 0.6 0.6],'MarkerSize',3)
           DVH= [DVH VO(t1,t2)];
        end        
        VPONTO=S'*PT*XT+XT'*PT*S+PHX;
        if VPONTO>0
           plot(x1(t1,t2),x2(t1,t2),'ko','MarkerSize',5)
           DVH= [DVH VO(t1,t2)];
        end            
     end
  end
LM=max(DVH)
LM=0.13; %0.13>LM anterior
% LM=0.161;
[C,h]=contour(x1,x2,VO,[LM,LM],'k-','LineWidth',2)
[C2,h2]=contour(x1,x2,VO,[0.49001,0.49001],'Color',[0.5 0.5 0.5],'LineWidth',2)

  [x1,x2] = meshgrid(dlim(1):.3:dlim(2),dlim(1):.3:dlim(2));
  [tF1,tF2]=size(x1);
  % Laço usado para criar o retrato de fase do sistema fuzzy
  VO=[];  DVH=[];
  LM=1.4*LM;
  for t1=1:tF1,
     for t2=1:tF2,
        XT=[x1(t1,t2);x2(t1,t2)]; % Vetor de estado.
        % Funções de pertinência
        h12=(x1(t1,t2)^2+x2(t1,t2)^2)/50; h11=1-h12;
        h21=(x2(t1,t2)^2)/25;  h22=1-h21;
        Ah1=h11*A{1,1}+h12*A{1,2};
        Ah2=h21*A{2,1}+h22*A{2,2};        
        S1=Ah1*XT; % Subsistema 1
        S2=Ah2*XT; % Subsistema 2         
%         dh12=(x1(t1,t2)*S1(1)+x2(t1,t2)*S1(2))/25; % h12ponto
%         dh21=2*x2(t1,t2)*S2(2)/25;
        dh12=(x1(t1,t2)*(alfa(1)*S1(1)+alfa(2)*S2(1))+x2(t1,t2)*(alfa(1)*S1(2)+alfa(2)*S2(2)))/25; % h12ponto
        dh21=2*x2(t1,t2)*(alfa(1)*S1(2)+alfa(2)*S2(2))/25;
        PT = h12*P{1,2}+h21*P{2,1};
        VO =  XT'*PT*XT;
        DPh=dh12*P{1,2}+dh21*P{2,1};
        DV1=XT'*DPh*XT+2*XT'*PT*S1; % Derivada da V no subsistema 1
        DV2=XT'*DPh*XT+2*XT'*PT*S2; % Derivada da V no subsistema 2
        if (VO > LM)  & (VO < 0.49001)
%            plot(x1(t1,t2),x2(t1,t2),'g*','MarkerSize',3);
%            plot(x1(t1,t2),x2(t1,t2),'*','Color',[0.6 0.6 0.6],'MarkerSize',3)
           if  DV1 < DV2
                dx1(t1,t2)=S1(1);
                dx2(t1,t2)=S1(2);
           else
                dx1(t1,t2)=S2(1);
                dx2(t1,t2)=S2(2);
           end
        else
                x1(t1,t2)=0;
                x2(t1,t2)=0;
                dx1(t1,t2)=0;
                dx2(t1,t2)=0;
        end
     end
  end

  sz = sqrt(dx1.^2 + dx2.^2); % The length of each arrow.
  dx1=dx1./sz;
  dx2=dx2./sz;
  quiver(x1,x2,dx1,dx2,1,'-','Color',[0.1 .8 0.2])
  xlabel('x1','fontsize',16);
  ylabel('x2','fontsize',16);
  text(0.6038,1.4054,'fi1','fontsize',14,'fontname','helvetica')
  text(1.8055,1.76,'fi2','fontsize',14,'fontname','helvetica')
  text(-0.3289,0.8889,'Ca','fontsize',14,'fontname','helvetica')

  
  % % CONDIÇÕES INICIAIS USADAS NAS SIMULAÇÕES DOS SUBSISTEMAS
X01=[-3;-1];
X02=[3.5;0.8];
%% Simulando o sistema fuzzy chaveado para a condição inicial X01
TF=5;
Tempo=[0:.001:TF];
Laux=[];
Taux=[];
CHMIN = 100; % Número máximo de iterações que um sistema pode ficar ativo
step=1;
global xsw at CHR1 CHR2;
at=1;  % Variavel diz qual o sistema q está ativo;
xsw=1; % Variavel controla o numero maximo de iterações que um sistema permanece ativo;
CHR1=CHMIN;
%rng(19);
CHR2=CHMIN;
save ex2.mat Laux Taux;
save ex2APMfi.mat A P LM CHMIN
y=ode1('exemplo_2_ODE',Tempo,X01);
St=length(Tempo);
figure(2);
plot(Tempo(1:step:St),y(1:step:St,1),'b-',...
     Tempo(1:step:St),y(1:step:St,2),'r--','LineWidth',2)
%% Gerando o gráfico da lei de chaveamento
load ex2.mat Laux Taux;
STaux=length(Taux);
hold on
plot(Taux(1:step:STaux),Laux(1:step:STaux),'-','Color',[0.5 0.5 0.5],'LineWidth',2)
%plot(Taux(1:step:STaux),Laux(1:step:STaux),'b','LineWidth',2)
xlabel('t','fontsize',16);
hleg= legend('X111(t)','X2','sig');
set(hleg,'fontsize',16,'fontname','helvetica');
axis([-.01 TF X01(1)-.2 4.2]);

figure(1)
 hold on
 plot(y(1:step:St,1),y(1:step:St,2),'b-','LineWidth',2)

 %% Gerando o gráfico da Função de Lyapunov
 for i=1:St,
    % Definindo as funções de pertinência do Subsistema 1
    h12=(y(i,1)^2+y(i,2)^2)/50;
    h11=1-h12;
    h21=(y(i,2)^2)/25;
    h22=1-h21;
    Vc(i)=[y(i,1) y(i,2)]*(h12*P{1,2}+h21*P{2,1})*[y(i,1);y(i,2)];
end
figure(3)
hold on
plot(Tempo,Vc,'b','LineWidth',2)
xlabel('t','fontsize',16)
ylabel('Vch','fontsize',16)
axis([-.02 TF 0 0.37])

%% Simulando o sistema fuzzy chaveado para a condição inicial X02
Laux=[]; Taux=[];
save ex2.mat Laux Taux;
y=ode1('exemplo_2_ODE',Tempo,X02);
St=length(Tempo);
figure(4);
plot(Tempo(1:step:St),y(1:step:St,1),'b-',...
     Tempo(1:step:St),y(1:step:St,2),'r--','LineWidth',2)
%% Gerando o gráfico da lei de chaveamento
load ex2.mat Laux Taux;
STaux=length(Taux);
hold on
% plot(Taux(1:step:STaux),Laux(1:step:STaux),'b','LineWidth',2)
plot(Taux(1:step:STaux),Laux(1:step:STaux),'-','Color',[0.5 0.5 0.5],'LineWidth',2)
xlabel('t','fontsize',16);
hleg= legend('X111(t)','X2','sig');
set(hleg,'fontsize',16,'fontname','helvetica');
axis([-.01 TF -.7 3.51]);

figure(1)
 hold on
 plot(y(1:step:St,1),y(1:step:St,2),'r-','LineWidth',2)

 %% Gerando o gráfico da Função de Lyapunov
 for i=1:St,
    % Definindo as funções de pertinência do Subsistema 1
    h12=(y(i,1)^2+y(i,2)^2)/50;
    h11=1-h12;
    h21=(y(i,2)^2)/25;
    h22=1-h21;
    Vc(i)=[y(i,1) y(i,2)]*(h12*P{1,2}+h21*P{2,1})*[y(i,1);y(i,2)];
end
figure(5)
hold on
plot(Tempo,Vc,'r','LineWidth',2)
xlabel('t','fontsize',16)
ylabel('Vch','fontsize',16)
axis([-.02 TF 0 0.45])

% Autovalor_subsitema1 = [eig(A{1,1}) eig(A{1,2})]
% Autovalor_subsitema2 = [eig(A{2,1}) eig(A{2,2})]
else
     display('Lmis infactiveis')
end