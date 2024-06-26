% To display a curve obtained using de Casteljau subdivision
% Used by run_decas-subdiv_g1 (or g2)
% returns the x and y coordinates of the sequence of points
% as row vectors
% The polyline is not plotted by this function
% n = level of recursion
% 
%  cpoly is a 2 x (m+1) matrix whose first row consists of x-coordinates
%  and second row of y-coordinates of m + 1 control points

%showd_decas_subdiv2 take as input a control polygon cpoly and output the x and y coordinates of the cubic
function [x, y] = show_decas_subdiv2(cpoly,n)
    r = 0; s= 1;
    lpoly = itersubdiv(cpoly,n);
    lnodes = makelist(lpoly);
    x = lnodes(1,:);
    y = lnodes(2,:);
end