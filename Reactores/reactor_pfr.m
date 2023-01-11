fprintf('Prueba 3er parcial sección Ejercicios\n')
fprintf('\nNombre: Luis Vargas')
Temp=308.15;
x=[0.4;0.6];
fprintf('Moléculas involucradas\nCH3-CH2NH-CH2-CH3 (1)\nCH3-(CH2)5-CH3 (2)\n\n')
fprintf('\n\nTabla con valores tomados de bibliografía\n')
T=table({'CH3'; 'CH2'; 'CH2NH'},[1;2;33],[0.9011;0.6744;1.2070],[0.848;0.540;0.936],[2;1;1],[2;5;0]);
T.Properties.VariableNames={'Moléculas','k','R_K','Q_K','v_k^(1)','v_k^(2)'}
fprintf('valores de r')
r1=T.("v_k^(1)")(1)*T.("R_K")(1)+T.("v_k^(1)")(2)*T.("R_K")(2)+T.("v_k^(1)")(3)*T.("R_K")(3);
r2=T.("v_k^(2)")(1)*T.("R_K")(1)+T.("v_k^(2)")(2)*T.("R_K")(2)+T.("v_k^(2)")(3)*T.("R_K")(3);
ri=[r1;r2]
fprintf('valores de q')
q1=T.("v_k^(1)")(1)*T.("Q_K")(1)+T.("v_k^(1)")(2)*T.("Q_K")(2)+T.("v_k^(1)")(3)*T.("Q_K")(3);
q2=T.("v_k^(2)")(1)*T.("Q_K")(1)+T.("v_k^(2)")(2)*T.("Q_K")(2)+T.("v_k^(2)")(3)*T.("Q_K")(3);
q=[q1;q2]
fprintf('valores de e_ki')
e_k1=[(T.("v_k^(1)")(1)*T.("Q_K")(1))/q1;(T.("v_k^(1)")(2)*T.("Q_K")(2))/q1;(T.("v_k^(1)")(3)*T.("Q_K")(3))/q1];
e_k2=[(T.("v_k^(2)")(1)*T.("Q_K")(1))/q2;(T.("v_k^(2)")(2)*T.("Q_K")(2))/q2;(T.("v_k^(2)")(3)*T.("Q_K")(3))/q2];
e_ik=table([1;2;33],e_k1,e_k2);
e_ik.Properties.VariableNames={'k','i=1','i=2'}
fprintf('valores de a_mk')
a_1k=[0;0;255.7];
a_2k=[0;0;255.7];
a_33k=[65.33;65.33;0];
a_mk=[a_1k a_2k a_33k]
fprintf('valores de tau_mk')
tau_mk=exp(-a_mk./Temp)
fprintf('valores de B_ik')
B_i1=[e_k1(1)*tau_mk(1,1)+e_k1(2)*tau_mk(1,2)+e_k1(3)*tau_mk(1,3);e_k2(1)*tau_mk(1,1)+e_k2(2)*tau_mk(1,2)+e_k2(3)*tau_mk(1,3)];
B_i2=[e_k1(1)*tau_mk(2,1)+e_k1(2)*tau_mk(2,2)+e_k1(3)*tau_mk(2,3);e_k2(1)*tau_mk(2,1)+e_k2(2)*tau_mk(2,2)+e_k2(3)*tau_mk(2,3)];
B_i33=[e_k1(1)*tau_mk(3,1)+e_k1(2)*tau_mk(3,2)+e_k1(3)*tau_mk(3,3);e_k2(1)*tau_mk(3,1)+e_k2(2)*tau_mk(3,2)+e_k2(3)*tau_mk(3,3)];
B_ik=table([1;2],B_i1,B_i2,B_i33);
B_ik.Properties.VariableNames={'i','k=1','k=2','k=3'}
fprintf('valores de theta_i')
theta_1=(x(1)*q(1)*e_ik.("i=1")(1)+x(2)*q(2)*e_ik.("i=2")(1))/(x(1)*q(1)+x(2)*q(2));
theta_2=(x(1)*q(1)*e_ik.("i=1")(2)+x(2)*q(2)*e_ik.("i=2")(2))/(x(1)*q(1)+x(2)*q(2));
theta_3=(x(1)*q(1)*e_ik.("i=1")(3)+x(2)*q(2)*e_ik.("i=2")(3))/(x(1)*q(1)+x(2)*q(2));
theta_i=[theta_1;theta_2;theta_3]
fprintf('valores de s_k')
s_1=theta_1.*tau_mk(1,1)+theta_2.*tau_mk(1,2)+theta_3.*tau_mk(1,3);
s_2=theta_1.*tau_mk(2,1)+theta_2.*tau_mk(2,2)+theta_3.*tau_mk(2,3);
s_3=theta_1.*tau_mk(3,1)+theta_2.*tau_mk(3,2)+theta_3.*tau_mk(3,3);
s_k=[s_1;s_2;s_3]
fprintf('valores de J_i')
J1=r1/(r1*x(1)+r2*x(2));
J2=r2/(r1*x(1)+r2*x(2));
J_i=[J1;J2];
fprintf('valores de L_i')
L1=q1/(q1*x(1)+q2*x(2))
L2=q2/(q1*x(1)+q2*x(2))
L_i=[L1;L2];
fprintf('valores de (lnyi)^C')
lny1C=1-J1+log(J1)-5*q1*(1-J1/L1+log(J1/L1));
lny2C=1-J2+log(J2)-5*q2*(1-J2/L2+log(J2/L2));
lnyi_C=[lny1C;lny2C]
fprintf('valores de (lnyi)^R')
lny1R=q1*(1-((theta_1*B_ik.("k=1")(1)/s_1-e_ik.("i=1")(1)*log(B_ik.("k=1")(1)/s_1)+(theta_2*B_ik.("k=2")(1)/s_2-e_ik.("i=1")(2)*log(B_ik.("k=2")(1)/s_2))+(theta_3*B_ik.("k=3")(1)/s_3-e_ik.("i=1")(3)*log(B_ik.("k=3")(1)/s_3)))));
lny2R=q2*(1-((theta_1*B_ik.("k=1")(2)/s_1-e_ik.("i=2")(1)*log(B_ik.("k=1")(2)/s_1)+(theta_2*B_ik.("k=2")(2)/s_2-e_ik.("i=2")(2)*log(B_ik.("k=2")(2)/s_2))+(theta_3*B_ik.("k=3")(2)/s_3-e_ik.("i=2")(3)*log(B_ik.("k=3")(2)/s_3)))));
lnyi_C=[lny1R;lny2R]
y1=exp(lny1C+lny1R);
y2=exp(lny2C+lny2R);
fprintf('valores de gamma yi')
y=[y1;y2]
