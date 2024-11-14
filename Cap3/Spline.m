%Spline: Calcula los coeficientes de los polinomios de interpolación de
% grado d (1, 2, 3) para el conjunto de n datos (x,y), 
% mediante el método spline.
function [Tabla] = Spline(x, y, d)
    n = length(x);
    A = zeros((d + 1) * (n - 1));
    b = zeros((d + 1) * (n - 1), 1);
    cua = x .^ 2;
    cub = x .^ 3;
    
    %% lineal
    if d == 1
        c = 1;
        h = 1;
        for i = 1:n-1
            A(i, c) = x(i);
            A(i, c + 1) = 1;
            b(i) = y(i);
            c = c + 2;
            h = h + 1;
        end
        
        c = 1;
        for i = 2:n
            A(h, c) = x(i);
            A(h, c + 1) = 1;
            b(h) = y(i);
            c = c + 2;
            h = h + 1;
        end
        
    %% Cuadratic
    elseif d == 2
        c = 1;
        h = 1;
        for i = 1:n-1
            A(i, c) = cua(i);
            A(i, c + 1) = x(i);
            A(i, c + 2) = 1;
            b(i) = y(i);
            c = c + 3;
            h = h + 1;
        end
        
        c = 1;
        for i = 2:n
            A(h, c) = cua(i);
            A(h, c + 1) = x(i);
            A(h, c + 2) = 1;
            b(h) = y(i);
            c = c + 3;
            h = h + 1;
        end
        
        c = 1;
        for i = 2:n-1
            A(h, c) = 2 * x(i);
            A(h, c + 1) = 1;
            A(h, c + 3) = -2 * x(i);
            A(h, c + 4) = -1;
            b(h) = 0;
            c = c + 4;
            h = h + 1;
        end
        
        A(h, 1) = 2;
        b(h) = 0;
        
    %% Cubic
    elseif d == 3
        c = 1;
        h = 1;
        for i = 1:n-1
            A(i, c) = cub(i);
            A(i, c + 1) = cua(i);
            A(i, c + 2) = x(i);
            A(i, c + 3) = 1;
            b(i) = y(i);
            c = c + 4;
            h = h + 1;
        end
        
        c = 1;
        for i = 2:n
            A(h, c) = cub(i);
            A(h, c + 1) = cua(i);
            A(h, c + 2) = x(i);
            A(h, c + 3) = 1;
            b(h) = y(i);
            c = c + 4;
            h = h + 1;
        end
        
        c = 1;
        for i = 2:n-1
            A(h, c) = 3 * cua(i);
            A(h, c + 1) = 2 * x(i);
            A(h, c + 2) = 1;
            A(h, c + 4) = -3 * cua(i);
            A(h, c + 5) = -2 * x(i);
            A(h, c + 6) = -1;
            b(h) = 0;
            c = c + 4;
            h = h + 1;
        end
        
        c = 1;
        for i = 2:n-1
            A(h, c) = 6 * x(i);
            A(h, c + 1) = 2;
            A(h, c + 4) = -6 * x(i);
            A(h, c + 5) = -2;
            b(h) = 0;
            c = c + 4;
            h = h + 1;
        end
        
        A(h, 1) = 6 * x(1);
        A(h, 2) = 2;
        b(h) = 0;
        h = h + 1;
        A(h, c) = 6 * x(end);
        A(h, c + 1) = 2;
        b(h) = 0;
    end

    val = inv(A) * b;
    Tabla = reshape(val, d + 1, n - 1)';
    
    %% Graficar los polinomios spline
    hold on;
    for i = 1:n-1
        % Obtener los coeficientes de cada polinomio
        coef = Tabla(i, :);
        
        % Crear puntos en el intervalo [x(i), x(i+1)]
        x_vals = linspace(x(i), x(i + 1), 100);
        
        % Evaluar el polinomio en x_vals
        if d == 1
            y_vals = coef(1) * x_vals + coef(2);
        elseif d == 2
            y_vals = coef(1) * x_vals.^2 + coef(2) * x_vals + coef(3);
        elseif d == 3
            y_vals = coef(1) * x_vals.^3 + coef(2) * x_vals.^2 + coef(3) * x_vals + coef(4);
        end
        
        % Graficar el tramo del spline
        plot(x_vals, y_vals, 'b', 'LineWidth', 1.5);
    end
    
    % Graficar los puntos originales
    plot(x, y, 'ro', 'MarkerFaceColor', 'r');

    xlabel('x');
    ylabel('y');
    title('Spline de Interpolación');
    grid on;
    hold off;
end
