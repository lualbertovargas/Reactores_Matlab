% Definir la longitud del reactor (L) y el paso de tiempo (dt)
L = 3;
dt = 0.1;

% Definir el número de pasos de tiempo (n)
n = L/dt;

% Inicializar el vector de concentraciones (c)
c = zeros(n+1, 2);

% Asignar las concentraciones iniciales
c(1, 1) = 1.8;
c(1, 2) = 0.4;

% Inicializar el vector de flujos (F)
F = zeros(n, 2);

% Inicializar el vector de tiempos (t)
t = zeros(n+1, 1);
k1=1;
k2=1;
k3=1;
% Bucle de simulación
for i = 1:n
    % Calcular el flujo de entrada al reactor
    F(i, 1) = 12;
    F(i, 2) = 14;
    
    % Calcular las concentraciones en el reactor
    c(i+1,1)=c(i,1)+(F(i,1)-k1*c(i,1)-k2*c(i,1)*c(i,2))*dt/L;
    c(i+1,2)=c(i,2)+(F(i,2)-k3*c(i,1)*c(i,2))*dt/L;
    
    % Actualizar el tiempo
    t(i+1) = t(i) + dt;
end

% Gráfica de las concentraciones en función del tiempo
plot(t, c(:, 1), 'b', t, c(:, 2), 'r');
xlabel('Tiempo (s)');
ylabel('Concentracion (s)');
