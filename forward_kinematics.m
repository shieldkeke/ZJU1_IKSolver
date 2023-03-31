syms t1 t2 t3 t4 t5 t6;
alpha = [0,-sym(pi)/2, 0, 0, sym(pi)/2, sym(pi)/2];
a = [0, 0, 0.185, 0.17, 0, 0];
d = [0.23, -0.054, 0, 0.077,0.077, 0.0855];
theta = [t1, -sym(pi)/2 + t2, t3, sym(pi)/2 + t4, sym(pi)/2 + t5, t6];
product = eye(4);
for i = 1:6
    m = transMat(alpha(i), a(i), d(i), theta(i));
    product = product * m;
end
simplify(product)
for i=1:4
    for j=1:4
        simplify(product(i,j))
    end
end

table = [0,0,0,0,0,0;
         sym(pi)/6, 0, sym(pi)/6, 0, sym(pi)/3, 0;
         sym(pi)/6, sym(pi)/6, sym(pi)/3, 0, sym(pi)/3, sym(pi)/6;
         sym(pi)/2, 0, sym(pi)/2, -sym(pi)/3, sym(pi)/3, sym(pi)/6;
         -sym(pi)/6, -sym(pi)/6, -sym(pi)/3, 0, sym(pi)/12, sym(pi)/2;
         sym(pi)/12, sym(pi)/12, sym(pi)/12, sym(pi)/12, sym(pi)/12, sym(pi)/12];

for i = 1:6
    for j = 1:6
        eval(['t',num2str(j),'= table(i,j);']);
    end
    trans = simplify(subs(product));
    trans_n = vpa(trans,4);
    eul = vpa(rotm2eul(double(trans_n(1:3,1:3)), 'XYZ') * 180 / pi,4)
end 