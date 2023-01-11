function dx=concaidap(W,x)
 
dx=zeros(2,1); 
k1=12;  
Cao=0.1;
vo=7.15;
a=50;


 dx(1)=(k1)*(Cao^2).*((1-x).^2);
 dx(2)=(1-vo*a);
%[W,x]=ode45(@caidap,[2 3 5],[0]);
plot(W,x,"--")
grid on
xlabel('t')
ylabel('q');
title('EDOs')
 
end