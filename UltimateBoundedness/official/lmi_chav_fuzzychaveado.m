function P = lmi_chav_fuzzychaveado (A,alfa,beta,xv,gama)
nA = size(A{1,1},2);  % Determina a dimensão da matriz A
alfa1=alfa(1); alfa2=alfa(2);

e1=[1;0]; e2=[0;1]; % Usados para limitar a curva de nível \Omega_1

P12= sdpvar(nA,nA,'symmetric');
P21= sdpvar(nA,nA,'symmetric');
P11= sdpvar(nA,nA,'symmetric');
P22= sdpvar(nA,nA,'symmetric');
  M= sdpvar(2*nA,2*nA,'symmetric');
R11= sdpvar(nA,nA,'full');
R12= sdpvar(nA,nA,'full');
R21= sdpvar(nA,nA,'full');
R22= sdpvar(nA,nA,'full');
L11= sdpvar(nA,nA,'full');
L12= sdpvar(nA,nA,'full');
L21= sdpvar(nA,nA,'full');
L22= sdpvar(nA,nA,'full');
 
% Criando as variaveis Garfo_{pk-ij}
a11=alfa1*(L11*A{1,1}+A{1,1}'*L11');
a21=P11-L11'+alfa1*(R11*A{1,1});
a22=-R11-R11';
Ga1111=[a11 a21';a21 a22]; % Garfo_{11-11}

a11=alfa2*(L22*A{2,2}+A{2,2}'*L22');
a21=P22-L22'+alfa2*(R22*A{2,2});
a22=-R22-R22';
Ga2222=[a11 a21';a21 a22]; % Garfo_{22-22}

Qalfa= Ga1111+Ga2222+M;

a11=alfa1*(L12*A{1,1}+A{1,1}'*L12');
a21=P12-L12'+alfa1*(R12*A{1,1});
a22=-R12-R12';
Ga1112=[a11 a21';a21 a22]; % Garfo_{11-12}
 
a11=alfa1*(L21*A{1,1}+A{1,1}'*L21');
a21=P21-L21'+alfa1*(R21*A{1,1});
a22=-R21-R21';
Ga1121=[a11 a21';a21 a22]; % Garfo_{11-21}

a11=alfa1*(L12*A{1,2}+A{1,2}'*L12');
a21=P12-L12'+alfa1*(R12*A{1,2});
a22=-R12-R12';
Ga1212=[a11 a21';a21 a22]; % Garfo_{12-12}
 
a11=alfa1*(L21*A{1,2}+A{1,2}'*L21');
a21=P21-L21'+alfa1*(R21*A{1,2});
a22=-R21-R21';
Ga1221=[a11 a21';a21 a22]; % Garfo_{12-21}

a11=alfa2*(L12*A{2,1}+A{2,1}'*L12');
a21=P12-L12'+alfa2*(R12*A{2,1});
a22=-R12-R12';
Ga2112=[a11 a21';a21 a22]; % Garfo_{21-12}
 
a11=alfa2*(L21*A{2,1}+A{2,1}'*L21');
a21=P21-L21'+alfa2*(R21*A{2,1});
a22=-R21-R21';
Ga2121=[a11 a21';a21 a22]; % Garfo_{21-21}

a11=alfa2*(L12*A{2,2}+A{2,2}'*L12');
a21=P12-L12'+alfa2*(R12*A{2,2});
a22=-R12-R12';
Ga2212=[a11 a21';a21 a22]; % Garfo_{22-12}
 
a11=alfa2*(L21*A{2,2}+A{2,2}'*L21');
a21=P21-L21'+alfa2*(R21*A{2,2});
a22=-R21-R21';
Ga2221=[a11 a21';a21 a22]; % Garfo_{22-21}

Restr = [];
% PT=[P12+P21 zeros(nA,nA); zeros(nA,nA) zeros(nA,nA)];
PT=[eye(nA) zeros(nA,nA); zeros(nA,nA) zeros(nA,nA)];
warning('off','YALMIP:nonstrict');
% LMI (29)
Restr = [Restr P12>=(e1*e1')/(xv^2) P12>=(e2*e2')/(xv^2) ...
               P21>=(e1*e1')/(xv^2) P21>=(e2*e2')/(xv^2) ...
               P12<=gama*eye(2) P21<=gama*eye(2)];
if beta==1
% LMI (21)
    Restr = [Restr Ga1212+Qalfa<=0];
% LMI (22)
    Restr = [Restr Ga1112+Qalfa<=0 Ga1121+Qalfa<=0];
% LMI (24)
    Restr = [Restr Ga2121-Qalfa<=0];
% LMI (25)
    Restr = [Restr Ga2212-Qalfa<=0 Ga2221-Qalfa<=0];    
% LMI (28)
    Restr = [Restr Ga1221+Ga2112<=0];
else
% LMI (21)
    Restr = [Restr Ga2121+Qalfa<=0];
% LMI (22)
    Restr = [Restr Ga2212+Qalfa<=0 Ga2221+Qalfa<0];
% LMI (24)
    Restr = [Restr Ga1212-Qalfa<=0];
% LMI (25)
    Restr = [Restr Ga1112-Qalfa<=0 Ga1121-Qalfa<=0];
% LMI (29)
    Restr = [Restr Ga2112+Ga1221<=0];
end

% Configurando o Solver.
opts=sdpsettings;
% opts.solver='lmilab';
 opts.solver='sedumi';
opts.verbose=0;
% Resolvendo as LMIs
sol = solvesdp(Restr,[],opts);
p=min(checkset(Restr));
if p > 0
    % Encontra o valor numérico das matrizes Pis
    P12=double(P12);
    P11=double(P11);
    P22=double(P22);
    P21=double(P21);
    R12=double(R12);
    R11=double(R11);
    R22=double(R22);
    R21=double(R21);
    L12=double(L12);
    L11=double(L11);
    L22=double(L22);
    L21=double(L21);
    %gama=double(gama);
    
    P=cell(2,2);
    P{1,1}=P11; P{2,1}=P21; P{2,2}=P22; P{1,2}=P12;
else
    display('LMIs infactíveis')
    P=[];
end
