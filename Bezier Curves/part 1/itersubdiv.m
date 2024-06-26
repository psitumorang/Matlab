function lpoly = itersubdiv(poly,n)
    lpoly = poly;
    for i=1:n
        lpoly = subdivstep(lpoly);
    end
end