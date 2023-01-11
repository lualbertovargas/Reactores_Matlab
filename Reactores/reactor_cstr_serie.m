%REACTORES EN SERIE CSTR (TRES)
% Valores de los parámetros del modelo
TIME = 0;
CA1 = 0.4;
CA2 = 0.2;
CA3 = 0.1;
CA0 = 1.8;
TAU = 2;
K = 0.5;
DELTA = 0.1;
TPRINT = 0;

% Imprimir encabezado
fprintf('TIME--------CA1--------CA2--------CA3\n');

% Bucle principal
while TIME <= 3
  % Evaluar derivadas
  CA1DOT = (CA0 - CA1) / TAU - K * CA1;
  CA2DOT = (CA1 - CA2) / TAU - K * CA2;
  CA3DOT = (CA2 - CA3) / TAU - K * CA3;

  % Imprimir resultados
  if TIME >= TPRINT
    fprintf('%8.3f %8.3f %8.3f %8.3f\n', TIME, CA1, CA2, CA3);
    TPRINT = TPRINT + 0.1;
    plot(TIME, CA1,'.')
    plot(TIME, CA2,'.')
    plot(TIME, CA3,'.')
    hold on 
    grid on   
    
  end

  % Actualizar variables
  CA1 = CA1 + CA1DOT * DELTA;
  CA2 = CA2 + CA2DOT * DELTA;
  CA3 = CA3 + CA3DOT * DELTA;
  TIME = TIME + DELTA;
end


