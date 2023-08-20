function yc = Example3_ODE (t, x)
% SCRIPT FOR ODE SOLVERS

load Example3_parameters.mat A P LM
load Example3.mat Laux Taux;
global ip MAX cont  % Global variables for switching law
  w01=(16-x(1)^2)/(149);     w11=1-w01;
  w02=(32-x(2)^2)/(198);     w12=1-w02;
h=zeros(1,4);
h(1)=w01*w02; h(2)=w01*w12; h(3)=w11*w02; h(4)=w11*w12;
% Calculating the matriz of each fuzzy subsystem
Ah1=h(1)*A{1,1}+h(2)*A{1,2}+h(3)*A{1,3}+h(4)*A{1,4};
Ah2=h(1)*A{2,1}+h(2)*A{2,2}+h(3)*A{2,3}+h(4)*A{2,4};
% Parametric matrix used in the Lyapunov function
 Ph=h(2)*(P{1,2}+P{2,2})+h(4)*(P{1,4}+P{2,4});
 Vh=x'*Ph*x;
% Derivative of V along of subsystem 1
VP1=x'*(Ah1'*Ph+Ph*Ah1)*x;

% Derivative of V along of subsystem 2
VP2=x'*(Ah2'*Ph+Ph*Ah2)*x;

% Checking if switching law (3) must be activated
if Vh>2*LM
    if VP1<0
       yc=Ah1*x;
       Laux=[Laux 1];
    else
       yc=Ah2*x;
       Laux=[Laux 2];
    end
else
% Any dwell time switching law can be used
switch ip
        case 1
            if(cont<=MAX)
                cont=cont+1;
                yc=Ah1*x;
                Laux=[Laux 1];
            else
                cont=1;
                ip=2;
                yc=Ah2*x;
                Laux=[Laux 2];
            end
        case 2
            if(cont<=MAX)
                cont=cont+1;
                yc=Ah2*x;
                Laux=[Laux 2];
            else
                cont=1;
                ip=1;
                yc=Ah1*x;
                Laux=[Laux 1];
            end
  end
end
Taux=[Taux t];
save Example3.mat Laux Taux;