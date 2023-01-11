function dx=sincaidapr(W,x)
 
dx=zeros(1,1); 
k1=12;  
Cao=0.1;
vo=7.15;

 dx(1)=k1*Cao*vo^(-1.*(1-x)^2);
%[W,x]=ode45(@caidap,[2 3 5],[0]);
plot(W,x,"--")
grid on
xlabel('t')
ylabel('q');
title('EDOs')
 
end