function [T] = newton(x0, Tol, niter)
    format long
    syms x

    f = (-exp(-x)) + x * (-1 + x) - x^(2/3)-9; % Definición de la función
    df = diff(f); % Derivada de la función
    c = 0; % Contador de iteraciones
    fm(c+1) = eval(subs(f, x0)); % Valor inicial de f(x0)
    fe = fm(c+1); % f en el punto x0
    dfm(c+1) = eval(subs(df, x0)); % Derivada f'(x0)
    dfe = dfm(c+1); % Almacena derivada en el punto x0
    E(c+1) = Tol + 1; % Error inicial
    error = E(c+1); % Almacena el valor del error
    xn(c+1) = x0; % Almacena el valor inicial de la solución aproximada
    
    % Inicializar matrices para la tabla
    n_vals = [];   % Para el número de iteraciones
    xn_vals = [];  % Para las aproximaciones
    fm_vals = [];  % Para los valores de la función
    dfm_vals = []; % Para los valores de la derivada
    E_vals = [];   % Para los errores
    
    % Almacenar la primera iteración
    n_vals = [n_vals; c];
    xn_vals = [xn_vals; x0];
    fm_vals = [fm_vals; fe];
    dfm_vals = [dfm_vals; dfe];
    E_vals = [E_vals; error];
    
    % Mostrar valores iniciales
    fprintf('Iteración %d: xn = %f, f(xn) = %f, f''(xn) = %f, Error = %f\n', c, x0, fe, dfe, error);
    
    while error > Tol && c < niter
        xn(c+2) = x0 - fe/dfe; % Siguiente aproximación
        fm(c+2) = eval(subs(f, xn(c+2))); % Valor de la función en xn
        fe = fm(c+2); % Actualiza el valor de f
        dfm(c+2) = eval(subs(df, xn(c+2))); % Derivada en xn
        dfe = dfm(c+2); % Actualiza derivada
        
        % Cálculo del error relativo
        if abs(xn(c+2)) > eps  % Evita división por cero
            E(c+2) = abs(xn(c+2) - x0) / abs(xn(c+2));
        else
            E(c+2) = abs(xn(c+2) - x0); % Fallback para evitar división por cero
        end
        
        error = E(c+2); % Actualiza el valor del error
        x0 = xn(c+2); % Actualiza xn
        c = c + 1; % Incrementa contador
        
        % Almacenar los valores de cada iteración
        n_vals = [n_vals; c];
        xn_vals = [xn_vals; x0];
        fm_vals = [fm_vals; fe];
        dfm_vals = [dfm_vals; dfe];
        E_vals = [E_vals; error];
        
        % Mostrar valores de cada iteración
        fprintf('Iteración %d: xn = %f, f(xn) = %f, f''(xn) = %f, Error = %f\n', c, x0, fe, dfe, error);
    end
    
    % Verifica si se cumple alguna condición de parada
    if fe == 0
        fprintf('%f es raíz de f(x) \n', x0);
    elseif error < Tol
        fprintf('%f es una aproximación de una raíz de f(x) con una tolerancia = %f \n', x0, Tol);
    elseif dfe == 0
        fprintf('%f es una posible raíz múltiple de f(x) \n', x0);
    else
        fprintf('Fracasó en %f iteraciones \n', niter);
    end
    
    % Crear la tabla con los datos solicitados
    T = table(n_vals, xn_vals, fm_vals, E_vals, ...
            'VariableNames', {'n', 'xn', 'fm', 'E'});
    
    % Mostrar la tabla final
    disp(T);
end
