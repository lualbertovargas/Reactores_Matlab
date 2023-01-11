clc
fprintf(2,'\t\t\t\tEncontrar erfc ')
syms x
y = input('\nIngrese el valor de f(kT/E(AB)): ');
fprintf('El valor de (kT/E(AB)) es:  ')
shortEcuation=(3.2995058+(2.1450892*x))/(1+(6.6586169*x)+(0.01425681*x^2));
format short
result2=solve(shortEcuation==y,x);

result2;
format short
res1=result2(1);
res2=result2(2);

x1=double(res1)
x2=double(res2)