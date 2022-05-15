function [F,properties,arguments]=max_model(X,method,options,extstruct)

switch method
    case 'graph'
        arguments=[];
        if length(extstruct.arg) == 1
            basis = getbase(extstruct.arg{1});
            inf_row = find(basis(:,1) == -inf);
            if length(inf_row)>0
                extstruct.arg{1}(inf_row) = [];
            end
            F = set(extstruct.var - extstruct.arg{1});
            arguments = extstruct.arg{1}(:);
        else
            arguments=[];
            F = set([]);
            for j = 1:length(extstruct.arg)
                basis = getbase(extstruct.arg{j});
                inf_row = find(basis(:,1) == -inf);
                if length(inf_row)>0
                    extstruct.arg{j}(inf_row) = [];
                end
                F = F + set(extstruct.var - extstruct.arg{j});
                arguments = [arguments;extstruct.arg{j}(:)];
            end
        end
        properties = struct('convexity','convex','monotonicity','increasing','definiteness','none');
    case 'milp'
        arguments = [];
        F = set([]);
        t = extstruct.var;
        for j = 1:length(extstruct.arg) % MAX(x,y)
            basis = getbase(extstruct.arg{j});
            inf_row = find(basis(:,1) == -inf);
            if length(inf_row)>0
                extstruct.arg{j}(inf_row) = [];
            end
            X = extstruct.arg{j};
            X = reshape(X,length(X),1);
            if prod(size(X)) == 1
                F = F + set(X == t);
            else
                [M,m] = derivebounds(X);
                %M = min(100,M);
                %m = max(-100,m);
                n = length(X);
                d = binvar(n,1);
                F = F + set(sum(d)==1);
                F = F + set(-(max(M)-min(m))*(1-d) <= t-X <= (max(M)-min(m))*(1-d));
                kk = [];
                ii = [];
                for i = 1:n
                    k = [1:1:i-1 i+1:1:n]';
                    ii = [ii;repmat(i,n-1,1)];
                    kk = [kk;k];
                    Mm = M(k)-m(i);
                end
                xii = extsubsref(X,ii);
                dii = extsubsref(d,ii);
                xkk = extsubsref(X,kk);
                F = F + set(xkk <= xii+(M(kk)-m(ii)).*(1-dii));
            end
            arguments = [arguments;extstruct.arg{j}(:)];
        end
        properties = struct('convexity','exact','monotonicity','exact','definiteness','none');

    otherwise
        F = [];
        return
end