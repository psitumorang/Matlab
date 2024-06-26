function [ud, ld] = subdecas(cpoly)

t = 0.5;

subdivs = cell(1,10);
subdivs{1} = cpoly;

for i = 1:size(cpoly,2)
    subdiv = [];    
    for j = 1:size(cpoly, 2)-1
        x_midpt = (1-t) * cpoly(1,j) + t * cpoly(1,j+1);
        y_midpt = (1-t) * cpoly(2,j) + t * cpoly(2,j+1);

        subdiv(1, j) = x_midpt;
        subdiv(2, j) = y_midpt;
    end

    subdivs{i+1} = subdiv;
    cpoly = subdiv;
end

subdivs = subdivs(~cellfun(@isempty, subdivs));
ud = [];
ld = [];

for i=1:size(subdivs,2)
    len = size(subdivs{i},2);
   
    ud = horzcat(ud,subdivs{i}(1:2,1));
    ld = horzcat(ld,subdivs{i}(1:2,len));
end

ld = fliplr(ld);
end

%disp(subdecas([0 1 2 3; 0 4 5 0]));