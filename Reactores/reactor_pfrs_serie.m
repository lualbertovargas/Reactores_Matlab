%la perturbación es la oportunidad de paso en la composicion del alimento
%en el tiempo cero de 0.5 a 0.55
%condiciones iniciales
CA=0.245;
T=600;
TJ=594.59;
V=48;
TIME=0;
VC=V*CA;
VT=V*T;
%VALORES DE PARAMETROS
TJ0=530;
F0=40;
T0=530;
CA0=0.5;
KC=4;
DELTA=0.0000001;
TPRINT=0;
fprintf('----TIME------CA-------T------V---------F------TJ----------FJ\n');
%Bucle principal
while TIME<=4.01
%DISTURBIO
CAO=0.55;
%CONTROLADORES DE RETROALIMENTACION
FJ=49.9-(KC*(600.-T));
F=40-(10.*(48-V));
%VELOCIDAD DE REACCION
K=(7.08E10)*exp(-30000./1.99*T);
Q=150*250*(T-TJ);
%EVALUAR TODO LOS DERIVADOS
VDOT=F0-F;
VCDOT=(F0*CAO)-(F*CA)-(V*K*CA);
VTDOT=(F0*T0)-(F*T)+(30000*V*K*CA)-(Q/(0.75*50));
TJDOT=((FJ*(TJ0-TJ))/3.85)+(Q/240);
%PARA IMPRIMIR LOS RESULTADOS
if TIME>=TPRINT
    fprintf('%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', TIME,CA,T,V,F,TJ,FJ);
    TPRINT=TPRINT+0.2;
    plot(TIME,CA,'r')
    plot(TIME, CA,'.')
    plot(TIME, T, '.')
    plot(TIME, V, '.')
    plot(TIME, F, '.')
    plot(TIME, TJ, '.')
    plot(TIME, FJ, '.')
    hold on
    grid on
end 
%Aplicando el delta a la variable
V=V+VDOT*DELTA;
VC=VC+VCDOT*DELTA;
VT=VT+VTDOT*DELTA;
TJ=TJ+TJDOT*DELTA;
TIME=TIME+DELTA;
CA=VC/V;
T=VT/V;
end