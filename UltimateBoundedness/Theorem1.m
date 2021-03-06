N = 2; %number of switching systems
n = 2; %state size
xv = [5,5]; %abs value of state limits in Z



calP = 1:N; %Set of the switching systems
calR = cell(1,length(calP)); %calR{p}, number of fuzzy subsystems of every switching system
calG = cell(1,length(calP)); %calG{p} subsets of calR{p} given

%% Create LMIs variables
for p = calP
    for k = calG{p}
        P{p,k} = sdpvar(nA,nA,'symmetric');
        L{p,k} = sdpvar(nA,nA,'full');
        R{p,k} = sdpvar(nA,nA,'full');
    end
end
M= sdpvar(2*nA,2*nA,'symmetric');

Upsilon = @(p,k,i,j)[alfa(i)*(L{p,k}*A{i,j} + A{i,j}'*L{p,k}'), (P{p,k}-L{p,k}'+alfa(i)*R{p,k}*A{i,j})';
                     P{p,k}-L{p,k}'+alfa(i)*R{p,k}*A{i,j},          -R{p,k}-R{p,k}'];

Q = 