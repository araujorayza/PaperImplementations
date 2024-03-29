function P = Example2_LMI_Theorem1 (A,alfa,xv,b,l,eps,b1,mu,be)
nA = size(A{1,1},2);  % Dimension of matrices A_pk
Ns = size(A,1);       % Number of subsystems
alfa1=alfa(1); alfa2=alfa(2);
e1=[1;0]; e2=[0;1];   % Boundary for level curve \Omega_1
P11= sdpvar(nA,nA,'symmetric');
P12= sdpvar(nA,nA,'symmetric');
P21= sdpvar(nA,nA,'symmetric');
P22= sdpvar(nA,nA,'symmetric');

P112= sdpvar(nA,nA,'symmetric');
P121= sdpvar(nA,nA,'symmetric');

  M= sdpvar(2*nA,2*nA,'symmetric');
R11= sdpvar(nA,nA,'full');
R12= sdpvar(nA,nA,'full');
R21= sdpvar(nA,nA,'full');
R22= sdpvar(nA,nA,'full');
L11= sdpvar(nA,nA,'full');
L12= sdpvar(nA,nA,'full');
L21= sdpvar(nA,nA,'full');
L22= sdpvar(nA,nA,'full');

% Creating variables \Upsilon_{pk-ij}
a11=alfa1*(L11*A{1,1}+A{1,1}'*L11')+((eps)/Ns)*P11-(eps*l/(Ns*(b^2)*b1))*eye(nA);
a21=(1/Ns)*P11-(1/Ns)*L11'+alfa1*(R11*A{1,1});
a22=(1/Ns)*(-R11-R11');
Ga1111=[a11 a21';a21 a22]; % \Upsilon_{11-11}

a11=alfa2*(L22*A{2,2}+A{2,2}'*L22')+((eps)/Ns)*P22-(eps*l/(Ns*(b^2)*b1))*eye(nA);
a21=(1/Ns)*P22-(1/Ns)*L22'+alfa2*(R22*A{2,2});
a22=(1/Ns)*(-R22-R22');
Ga2222=[a11 a21';a21 a22]; % \Upsilon_{22-22}

Qalfa= Ga1111+Ga2222+M;

a11=alfa1*(L12*A{1,1}+A{1,1}'*L12')+((eps)/Ns)*P12-(eps*l/(Ns*(b^2)*b1))*eye(nA);
a21=(1/Ns)*P12-(1/Ns)*L12'+alfa1*(R12*A{1,1});
a22=(1/Ns)*(-R12-R12');
Ga1112=[a11 a21';a21 a22]; % \Upsilon_{11-12}
 
a11=alfa1*(L21*A{1,1}+A{1,1}'*L21')+((eps)/Ns)*P21-(eps*l/(Ns*(b^2)*b1))*eye(nA);
a21=(1/Ns)*P21-(1/Ns)*L21'+alfa1*(R21*A{1,1});
a22=(1/Ns)*(-R21-R21');
Ga1121=[a11 a21';a21 a22]; % \Upsilon_{11-21}

a11=alfa1*(L12*A{1,2}+A{1,2}'*L12')+((eps)/Ns)*P12-(eps*l/(Ns*(b^2)*b1))*eye(nA);
a21=(1/Ns)*P12-(1/Ns)*L12'+alfa1*(R12*A{1,2});
a22=(1/Ns)*(-R12-R12');
Ga1212=[a11 a21';a21 a22]; % \Upsilon_{12-12}
 
a11=alfa1*(L21*A{1,2}+A{1,2}'*L21')+((eps)/Ns)*P21-(eps*l/(Ns*(b^2)*b1))*eye(nA);
a21=(1/Ns)*P21-(1/Ns)*L21'+alfa1*(R21*A{1,2});
a22=(1/Ns)*(-R21-R21');
Ga1221=[a11 a21';a21 a22]; % \Upsilon_{12-21}

a11=alfa2*(L12*A{2,1}+A{2,1}'*L12')+((eps)/Ns)*P12-(eps*l/(Ns*(b^2)*b1))*eye(nA);
a21=(1/Ns)*P12-(1/Ns)*L12'+alfa2*(R12*A{2,1});
a22=(1/Ns)*(-R12-R12');
Ga2112=[a11 a21';a21 a22]; % \Upsilon_{21-12}
 
a11=alfa2*(L21*A{2,1}+A{2,1}'*L21')+((eps)/Ns)*P21-(eps*l/(Ns*(b^2)*b1))*eye(nA);
a21=(1/Ns)*P21-(1/Ns)*L21'+alfa2*(R21*A{2,1});
a22=(1/Ns)*(-R21-R21');
Ga2121=[a11 a21';a21 a22]; % \Upsilon_{21-21}

a11=alfa2*(L12*A{2,2}+A{2,2}'*L12')+((eps)/Ns)*P12-(eps*l/(Ns*(b^2)*b1))*eye(nA);
a21=(1/Ns)*P12-(1/Ns)*L12'+alfa2*(R12*A{2,2});
a22=(1/Ns)*(-R12-R12');
Ga2212=[a11 a21';a21 a22]; % \Upsilon_{22-12}
 
a11=alfa2*(L21*A{2,2}+A{2,2}'*L21')+((eps)/Ns)*P21-(eps*l/(Ns*(b^2)*b1))*eye(nA);
a21=(1/Ns)*P21-(1/Ns)*L21'+alfa2*(R21*A{2,2});
a22=(1/Ns)*(-R21-R21');
Ga2221=[a11 a21';a21 a22]; % \Upsilon_{22-21}

Restr = [];
% LMI (19)
  Restr = [P12<=mu*eye(2) P21<=mu*eye(2)];
% LMI (20)
  Restr = [Restr P12>=(e1*e1')/(xv^2) P12>=(e2*e2')/(xv^2)...
                 P21>=(e1*e1')/(xv^2) P21>=(e2*e2')/(xv^2)];
% LMI (21)
  Restr = [Restr P112-P12>=0  P121-P21>=0];
if be==1
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
% Setting LMI solver
opts=sdpsettings;
% opts.solver='lmilab';
 opts.solver='sedumi';
opts.verbose=0;
% Solving LMIs
sol = solvesdp(Restr,[],opts);
p=min(checkset(Restr));
if p > 0
    % Obtaining numerical values of matrices P_pk
    P{1,1}=double(P11);
    P{1,2}=double(P12);
    P{2,1}=double(P21);
    P{2,2}=double(P22);
    R12=double(R12);
    R11=double(R11);
    R22=double(R22);
    R21=double(R21);
    L12=double(L12);
    L11=double(L11);
    L22=double(L22);
    L21=double(L21);
else
    display('Unfeasible LMIs')
    P=[];
end