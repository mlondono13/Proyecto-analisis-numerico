import sympy as sp

def secante(x0, x1, Tol, niter):
    x = sp.symbols('x')
    f = 2 * x**3 - 2 * x - 5  # Definición de la función
    fm = [f.subs(x, x0).evalf(), f.subs(x, x1).evalf()]  # Valores iniciales de f(x0) y f(x1)
    E = [Tol + 1]  # Error inicial
    error = E[0]
    
    # Inicializar listas para la tabla
    n_vals = [0, 1]       # Número de iteraciones
    xn_vals = [x0, x1]    # Aproximaciones
    fm_vals = [fm[0], fm[1]]  # Valores de la función
    E_vals = [error, error]  # Errores
    
    # Mostrar la cabecera de la tabla
    print('n\txn\tfm\tE')
    
    # Mostrar las primeras dos iteraciones
    print(f'{0}\t{x0:.6f}\t{fm[0]:.6f}\t{error:.6f}')
    print(f'{1}\t{x1:.6f}\t{fm[1]:.6f}\t{error:.6f}')
    
    c = 1  # Contador de iteraciones

    while error > Tol and c < niter:
        # Aplicar el método de la secante
        xn = x1 - fm[1] * (x1 - x0) / (fm[1] - fm[0])
        fm.append(f.subs(x, xn).evalf())  # Valor de la función en xn
        error = abs(xn - x1)  # Calcula el error

        # Almacenar los valores de cada iteración
        c += 1
        n_vals.append(c)
        xn_vals.append(xn)
        fm_vals.append(fm[2])
        E_vals.append(error)
        
        # Mostrar valores de cada iteración
        print(f'{c}\t{xn:.6f}\t{fm[2]:.6f}\t{error:.6f}')
        
        # Actualiza x0 y x1 para la siguiente iteración
        x0, x1 = x1, xn
        fm = [fm[1], fm[2]]  # Actualiza los valores de f(x0) y f(x1)

    # Verificar condiciones de parada
    if fm[2] == 0:
        print(f'{xn} es raíz de f(x)')
    elif error < Tol:
        print(f'{xn} es una aproximación de una raíz de f(x) con una tolerancia = {Tol}')
    else:
        print(f'Fracasó en {niter} iteraciones')

    # Crear un diccionario para la "tabla" de resultados
    T = {
        'n': n_vals,
        'xn': xn_vals,
        'fm': fm_vals,
        'E': E_vals
    }
    
    return T
