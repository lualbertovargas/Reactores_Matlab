clear
clc
CA0=1.8;
%flujo vol asumido
v=1;
k=0.5;
%calculo de flujo
FA0=CA0/v;
%variación de concentración
for i=1:CA0*10
    CA(19-i)=i/10;
    x(19-i,:)=(CA0-CA(19-i))/CA0;
end
%calculo constante velocidad
rA=k.*CA;
%flujo
F_R=FA0./rA';
hold on
plot(x,F_R)
xlabel('x')
ylabel('FA0/rA')
x1=[0.77 0.89 0.94];
y=FA0./(k.*[0.4 0.2 0.1]);
plot(x1,y,'g*')
figure
hold off
T=table(x,F_R)
fprintf('Para tres CSTR')
hold on
plot(x,F_R)
xlabel('x')
ylabel('FA0/rA')
title('Volumen para CSTR')
y=0:1:9;
x0=0.77.*ones(10,1);
y1=0:0.385:0.77;
x1=9.*ones(3,1);
plot(x0,y,'r-',y1,x1,'r-')
y=0:1:18;
x0=0.89.*ones(19,1);
y2=0.77:0.06:0.89;
x2=18.*ones(3,1);
plot(x0,y,'g-',y2,x2,'g-',0.77*ones(10,1),9:1:18,'g-')
y=0:1:36;
x0=0.94.*ones(37,1);
y3=0.89:0.025:0.94;
x3=36.*ones(3,1);
plot(x0,y,'k-',y3,x3,'k-',0.89*ones(19,1),18:1:36,'k-')
hold off

%cálculo de volumen de CSTR
V_CSTR(1)=FA0/(k*0.4)*(0.77);
V_CSTR(2)=FA0/(k*0.2)*(0.89-0.77);
V_CSTR(3)=FA0/(k*0.1)*(0.94-0.89);
VT_CSTR=V_CSTR(1)+V_CSTR(2)+V_CSTR(3)
%Para tres PFR
figure
fprintf('Para tres PFR')
hold on
plot(x,F_R)
xlabel('x')
ylabel('FA0/rA')
title('Volumen para PFR')

%límite del area de cada reactor
y=0:1:9;
x=0.77.*ones(10,1);
plot(x,y,'r-')
y=0:1:18;
x=0.89.*ones(19,1);
plot(x,y,'g-')
y=0:1:36;
x=0.94.*ones(37,1);
plot(x,y,'k-')
hold off

%Cálculo de volumen de reactor PFR
V_PFR(1)=0.1111/3*(F_R(1)+4*F_R(8)+F_R(15));
V_PFR(2)=0.0556/3*(F_R(15)+4*F_R(16)+F_R(17));
V_PFR(3)=0.0278/3*(F_R(17)+4*27.0162+F_R(18));
VT_PFR=V_PFR(1)+V_PFR(2)+V_PFR(3)