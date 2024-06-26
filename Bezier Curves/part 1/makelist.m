function lnodes = makelist(lpoly)
    lnodes = [];

    for i=1:size(lpoly,3)
        if (i == 1)
            lnodes = [lnodes, lpoly(:,:,i)];
        end
        if (i > 1)
            a = lpoly(:,:,i);
            a(:, 1) = [];
            lnodes = [lnodes, a];
        end
    end
end