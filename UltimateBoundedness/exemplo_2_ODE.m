% SCRIPT USADO PELO ODE45 PARA SIMULAR O SISTEMA FUZZY
function yc = exemplo_2_ODE (t, x)
load ex2APMfi.mat A P LM CHMIN
global xsw at CHR1 CHR2;
% Definindo as funções de pertinência do Exemplo 1
h(1,2)=(x(1)^2+x(2)^2)/50;
h(1,1)=1-h(1,2);
h(2,1)=(x(2)^2)/25;
h(2,2)=1-h(2,1);
Ah1=h(1,1)*A{1,1}+h(1,2)*A{1,2};
Ah2=h(2,1)*A{2,1}+h(2,2)*A{2,2};
S1=Ah1*x; % Subsistema 1
S2=Ah2*x; % Subsistema 2
dh(1,2)=(x(1)*S1(1)+x(2)*S1(2))/25; % h12ponto
dh(2,1)=2*x(2)*S2(2)/25;  % h21ponto
%% % Determina se o sistema permanece ativo no subsistema 1
Ph=h(1,2)*P{1,2}+h(2,1)*P{2,1};
DPh=dh(1,2)*P{1,2}+dh(2,1)*P{2,1};
Vch=x'*Ph*x;
DV1=x'*DPh*x+2*x'*Ph*S1; % Derivada da V no subsistema 1
DV2=x'*DPh*x+2*x'*Ph*S2; % Derivada da V no subsistema 2
if Vch > LM
    xsw=1;
    if  DV1 < DV2
        at=1;
        yc = SwitchingLaw(t,S1,S2);
    else
        at=2;
        yc = SwitchingLaw(t,S1,S2);
    end
else    
	switch at        
        case 1
             if  xsw <= CHR2
                   yc = SwitchingLaw(t,S1,S2);
                   xsw = xsw + 1;
             else
                   at=2;
                   yc = SwitchingLaw(t,S1,S2);
                   xsw = 1;
                   aux=CHR2;
                   while(abs(CHR2-CHR1)<30) 
                       CHR2=randi(CHMIN);
                   end
                   CHR1=aux;
             end
        case 2
             if  xsw <= CHR2
                   yc = SwitchingLaw(t,S1,S2);
                   xsw = xsw + 1;
             else
                   at=1;
                   yc = SwitchingLaw(t,S1,S2);
                   xsw = 1;
                   aux=CHR2;
                   while(abs(CHR2-CHR1)<30) 
                       CHR2=randi(CHMIN);
                   end
                   CHR1=aux;
             end
    end
end

function yc = SwitchingLaw(t,S1,S2)
global at;
load ex2.mat Laux Taux;
switch at
   case 1
       yc = S1;
       Laux = [Laux;1];
   case 2
       yc = S2;
       Laux = [Laux;2];
   case 3
       % Erro na lei de chaveamento
       display('Existe t para o qual o sistema não atende a condição de chaveamento')
       B2; % Essa linha pausa a execução do MATLAB       
end
Taux = [Taux;t];
save ex2.mat Laux Taux;
