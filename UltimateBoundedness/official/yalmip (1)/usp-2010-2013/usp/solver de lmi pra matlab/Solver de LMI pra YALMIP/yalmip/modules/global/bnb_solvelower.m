function  output = bnb_solvelower(lowersolver,relaxed_p,upper,lower)

if all(relaxed_p.lb==relaxed_p.ub)
    x = relaxed_p.lb;
    res = resids(relaxed_p,x);
    if all(res>-1e-7)
        output.problem = 0;
    else
        output.problem = 1;
    end
    output.Primal = x;
    return
end

p = relaxed_p;
removethese = p.lb==p.ub;
if nnz(removethese)>0 & all(p.variabletype == 0) & isempty(p.evalMap)% ~isequal(lowersolver,'callfmincongp') & ~isequal(lowersolver,'callgpposy')

    
    if ~isinf(upper) & nnz(p.Q)==0
        p.F_struc = [p.F_struc(1:p.K.f,:);upper-p.f -p.c';p.F_struc(1+p.K.f:end,:)];
        p.K.l=p.K.l+1;
%     elseif ~isinf(upper) & eig(full(p.Q))>=0 
%         x = relaxed_p.x_min;
%         grad = 2*p.Q*x + p.c;
%         f0 = p.f + p.c'*x + x'*p.Q*x;
%         p.F_struc = [p.F_struc(1:p.K.f,:);upper-f0 -grad';p.F_struc(1+p.K.f:end,:)];
%         p.K.l=p.K.l+1;
    end
       
    
    ptemp = p;
    if ~isempty(p.F_struc)
        if ~isequal(p.K.l,0) & p.options.bnb.ineq2eq
            eliminated = p.F_struc(:,1+find(removethese));
            affected = find(any(p.F_struc(:,1+find(removethese)),2));
        end
        p.F_struc(:,1)=p.F_struc(:,1)+p.F_struc(:,1+find(removethese))*p.lb(removethese);
        p.F_struc(:,1+find(removethese))=[];

        
%         if p.K.f > 0
%             lb = p.lb;lb(removethese) = [];
%             ub = p.ub;ub(removethese) = [];
%             dummy = p.F_struc(1:p.K.f,2:end);
%             Candbounds = sum(dummy | dummy,2)==1;
%             [ii,jj,kk] = find(dummy(Candbounds,:))
%             derivedBounds = full(p.F_struc(Candbounds,1)./(-kk));
%             if any(derivedBounds < lb(jj)) | any(derivedBounds > ub(jj))
%                 output.Problem = 1;
%                 output.Primal = zeros(length(p.c),1);
%             else
%                 keepthese = setdiff(1:length(p.c),removethese);
%                 p.lb(keepthese(jj)) = derivedBounds;
%                 p.ub(keepthese(jj)) = derivedBounds;
%                 removethese = p.lb==p.ub;
%                 p = ptemp;
%                 if ~isequal(p.K.l,0) & p.options.bnb.ineq2eq
%                     eliminated = p.F_struc(:,1+find(removethese));
%                     affected = find(any(p.F_struc(:,1+find(removethese)),2));
%                 end
%                 p.F_struc(:,1)=p.F_struc(:,1)+p.F_struc(:,1+find(removethese))*p.lb(removethese);
%                 p.F_struc(:,1+find(removethese))=[];
% 
%             end
%         end
    end

    
    p.c(removethese)=[];
    if nnz(p.Q)>0
        p.c = p.c + 2*p.Q(find(~removethese),find(removethese))*p.lb(removethese);
        p.Q(:,find(removethese))=[];
        p.Q(find(removethese),:)=[];
    else
        p.Q = spalloc(length(p.c),length(p.c),0);
    end
    p.lb(removethese)=[];
    p.ub(removethese)=[];
    p.x0(removethese)=[];
    p.monomtable(:,find(removethese))=[];
    p.monomtable(find(removethese),:)=[];
    p.variabletype = []; % Reset, to messy to recompute
    if ~isequal(p.K.l,0) & p.options.bnb.ineq2eq
        Beq = p.F_struc(1:p.K.f,1);
        Aeq = -p.F_struc(1:p.K.f,2:end);
        B   = p.F_struc(1+p.K.f:p.K.l+p.K.f,1);
        A   = -p.F_struc(1+p.K.f:p.K.l+p.K.f,2:end);        
        affected = affected(affected <= p.K.f + p.K.l);
        affected = affected(affected > p.K.f) - p.K.f;
        aaa = zeros(p.K.l,1);aaa(affected) = 1;
        A1 = A(find(~aaa),:);
        B1 = B(find(~aaa),:);
        [A,B,Aeq2,Beq2,index] = ineq2eq(A(affected,:),B(affected));
        if ~isempty(index)
            p.F_struc = [Beq -Aeq;Beq2 -Aeq2;B1 -A1;B -A;p.F_struc(1+p.K.f + p.K.l:end,:)];
            p.K.f = length(Beq) + length(Beq2);
            p.K.l = length(B) + length(B1);
        end
    end

%      % for safety, remove trailing garbage
%      zerorows = find(sum(abs(p.F_struc(1:(p.K.f+p.K.l),:)),2) == 0)
%      if ~isempty(zerorows)
%          p.F_struc(zerorows,:) = [];
%          p.K.f = p.K.f - nnz(zerorows <= p.K.f);
%          p.K.l = p.K.l - nnz(zerorows > p.K.f);
%      end
         
%     p.lb(isinf(p.lb)) = -1000;
%     p.ub(isinf(p.ub)) = 1000;
    output = feval(lowersolver,p);
    x=relaxed_p.c*0;
    x(removethese)=relaxed_p.lb(removethese);
    x(~removethese)=output.Primal;
    output.Primal=x;           
else
    output = feval(lowersolver,p);
end



function [A, B, Aeq, Beq, ind_eq] = ineq2eq(A, B)
% Copyright is with the following author(s):
%
% (C) 2006 Johan Loefberg, Automatic Control Laboratory, ETH Zurich,
%          loefberg@control.ee.ethz.ch
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

[ne, nx] = size(A);
Aeq = [];
Beq = [];
ind_eq = [];
sumM = sum(A, 2) + B;
for ii = 1:ne-1,   
    s = sumM(1);
    
    % get matrix which contains all rows starting from ii+1 
    sumM = sumM(2:end,:);
    
    % possible candidates are those rows whose sum is equal to the sum of the
    % original row
    possible_eq = find(abs(sumM + s) < 1e-12);
    if isempty(possible_eq),
        continue
    end
    possible_eq = possible_eq + ii;
    b1 = B(ii);    
    a1 = A(ii,  :) ;
    
    % now compare if the two inequalities (the second one with opposite
    % sign) are really equal (hence they form an equality constraint)
    for jj = possible_eq',
        % first compare the B part, this is very cheap
        if abs(b1 + B(jj)) < 1e-12,
            % now compare the A parts as well
            if norm(a1 + A(jj,  :) , Inf) < 1e-12,
                % jj-th inequality together with ii-th inequality forms an equality
                % constraint
                ind_eq = [ind_eq; ii jj];
                break
            end
        end
    end
end

if isempty(ind_eq),
    % no equality constraints
    return
else
    % indices of remaining constraints which are inequalities
    ind_ineq = setdiff(1:ne, ind_eq(:));
    Aeq = A(ind_eq(:,1),  :) ;
    Beq = B(ind_eq(:,1));
    A = A(ind_ineq,  :) ;
    B = B(ind_ineq);
end