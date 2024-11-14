function [T] = secante(x0, x1, Tol, niter)
    syms x

    f = 2*x^3 - 2*x - 5; % Definición de la función
    fm(1) = eval(subs(f, x0)); % Valor inicial f(x0)
    fm(2) = eval(subs(f, x1)); % Valor inicial f(x1)
    E(1) = Tol + 1; % Error inicial
    error = E(1); % Almacena el valor del error
    
    % Inicializar matrices para la tabla
    n_vals = [];   % Para el número de iteraciones
    xn_vals = [];  % Para las aproximaciones
    fm_vals = [];  % Para los valores de la función
    E_vals = [];   % Para los errores
    
    % Mostrar la cabecera de la tabla
    fprintf('n,xn,fm,E\n');

    % Almacenar la primera iteración
    n_vals = [n_vals; 0];
    xn_vals = [xn_vals; x0];
    fm_vals = [fm_vals; fm(1)];
    E_vals = [E_vals; error];
    fprintf('%d,%.6f,%.6f,%.6f\n', 0, x0, fm(1), error);
    
    % Almacenar la segunda iteración
    n_vals = [n_vals; 1];
    xn_vals = [xn_vals; x1];
    fm_vals = [fm_vals; fm(2)];
    E_vals = [E_vals; error];
    fprintf('%d,%.6f,%.6f,%.6f\n', 1, x1, fm(2), error);

    c = 1; % Contador de iteraciones

    while error > Tol && c < niter
        % Aplicar el método de la secante
        xn = x1 - fm(2) * (x1 - x0) / (fm(2) - fm(1));
        fm(3) = eval(subs(f, xn)); % Valor de la función en xn
        error = abs(xn - x1); % Calcula el error
        
        % Almacenar los valores de cada iteración
        c = c + 1; % Incrementa el contador
        n_vals = [n_vals; c];
        xn_vals = [xn_vals; xn];
        fm_vals = [fm_vals; fm(3)];
        E_vals = [E_vals; error];
        
        % Mostrar valores de cada iteración
        fprintf('%d,%.6f,%.6f,%.6f\n', c, xn, fm(3), error);
        
        % Actualiza x0 y x1 para la siguiente iteración
        x0 = x1;
        x1 = xn;
    end
    
    % Verifica si se cumple alguna condición de parada
    if fm(3) == 0
        fprintf('%f es raíz de f(x) \n', xn);
    elseif error < Tol
        fprintf('%f es una aproximación de una raíz de f(x) con una tolerancia = %f \n', xn, Tol);
    else
        fprintf('Fracasó en %d iteraciones \n', niter);
    end
    
    % Crear la tabla con los datos solicitados
    T = table(n_vals, xn_vals, fm_vals, E_vals, ...
              'VariableNames', {'n', 'xn', 'fm', 'E'});
end
