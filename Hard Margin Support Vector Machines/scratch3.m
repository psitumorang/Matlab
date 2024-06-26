v = cell(2, 50);
u = cell(2, 50);

rng(14175332);

u{4} = 10.1*randn(2,50)+18;
v{4} = -10.1*randn(2,50)-18;

u = u{4};
v = v{4};

rho = 10;

p = size(u,2); q = size(v,2); n = size(u,1);
[A,c,X,Pa,qa] = buildhardSVM2(u,v);

tolr = 10^(-10); tols = 10^(-10); iternum = 180000;
[lam,U,nr,ns,kk] = qsolve1(Pa, qa, A, c, rho, tolr, tols, iternum);
if kk > iternum
   fprintf('** qsolve did not converge. Problem not solvable ** \n')
end
lamb = lam(1:p,1); 
mu = lam(p+1:p+q,1);

% calculate w

ui = [];
sum_i = [0; 0];
for i = 1:p
    ui_x = lamb(i) * u(1, i);
    ui_y = lamb(i) * u(2, i);

    ui(1,i) = ui_x;
    ui(2,i) = ui_y;
    
    sum_i(1) = sum_i(1) + ui_x;
    sum_i(2) = sum_i(2) + ui_y;

end

vj = [];
sum_j = [0; 0];

for j = 1:q
    vj_x = mu(j) * v(1, j);
    vj_y = mu(j) * v(2, j);

    vj(1,j) = vj_x;
    vj(2,j) = vj_y;
    
    sum_j(1) = sum_j(1) + vj_x;
    sum_j(2) = sum_j(2) + vj_y;
end


w = sum_i - sum_j;

% calculate numsvl1, numsvm1

numsvl1 = 0;
numsvm1 = 0;

tol = 10 ^(-10);

for i = 1:p
    if lamb(i) > tol
        numsvl1 = numsvl1 +1;
    end
end

for j = 1:q
    if mu(j) > tol
        numsvm1 = numsvm1 +1;
    end
end

% calculate b

sum_ui = [0; 0];
sum_vj = [0; 0];

for i = 1:p
    if lamb(i) > 0
        sum_ui(1) = sum_ui(1) + u(1,i);
        sum_ui(2) = sum_ui(2) + u(2,i);
    end
end

for j = 1:q
    if mu(j) > 0
        sum_vj(1) = sum_vj(1) + v(1,j);
        sum_vj(2) = sum_vj(2) + v(2,j);
    end
end

b = transpose(w) * ((sum_ui)/numsvl1 + (sum_vj)/numsvm1) /2;
