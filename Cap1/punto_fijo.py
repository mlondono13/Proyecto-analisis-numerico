import sympy as sp
import pandas as pd

def punto_fijo(x0, Tol, niter, A):
    x = sp.symbols('x')
    # Definir la función y su iteración g(x) para el método de punto fijo
    f = sp.sin(x - A * 10**(-3)) - x
    g = sp.sin(x - A * 10**(-3))

    # Inicializar variables
    c = 0
    fm = [float(f.subs(x, x0))]  # Evaluación inicial de f(x0)
    fe = fm[c]  # Valor de f(x0)
    E = [Tol + 1]  # Error inicial
    error = E[c]
    xn = [x0]  # Almacena las aproximaciones de la raíz
    N = [c]  # Contador de iteraciones

    # Ciclo de iteración de punto fijo
    while error > Tol and fe != 0 and c < niter:
        xn.append(float(g.subs(x, x0)))  # Siguiente aproximación
        fm.append(float(f.subs(x, xn[c + 1])))  # Evaluación de f en xn
        fe = fm[c + 1]  # Actualizar el valor de f
        E.append(abs(xn[c + 1] - x0))  # Cálculo del error absoluto
        error = E[c + 1]  # Actualizar error
        x0 = xn[c + 1]  # Actualizar xn
        N.append(c + 1)  # Incrementar iteración
        c += 1

    # Verificar condición de parada
    if fe == 0:
        print(f"{x0} es raíz de f(x)")
    elif error < Tol:
        print(f"{x0} es una aproximación de una raíz de f(x) con una tolerancia = {Tol}")
    else:
        print(f"Fracasó en {niter} iteraciones")

    # Crear la tabla de resultados con pandas
    T = pd.DataFrame({
        'n': N,
        'xn': xn,
        'fm': fm,
        'E': E
    })

    # Mostrar la tabla final
    print(T)
    return T

# Ejemplo de uso
T = punto_fijo(1.0, 1e-5, 50, 2)
