function display(X)
%display           Overloaded

% Author Johan L�fberg 
% $Id: display.m,v 1.1 2006/08/24 17:12:27 joloef Exp $  

disp(['Optimization object with ' num2str(X.n) ' inputs and ' num2str(length(X.map)) ' outputs.'])
