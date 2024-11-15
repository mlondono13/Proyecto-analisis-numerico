import sympy as sp
import pandas as pd

def newton(x0, Tol, niter):
    x = sp.symbols('x')
    # Definir la función y su derivada
    f = -sp.exp(-x) + x * (-1 + x) - x**(2/3) - 9
    df = sp.diff(f, x)  # Derivada de la función

    # Inicializar variables
    c = 0
    fm = [float(f.subs(x, x0))]
    fe = fm[c]
    dfm = [float(df.subs(x, x0))]
    dfe = dfm[c]
    E = [Tol + 1]  # Error inicial
    error = E[c]
    xn = [x0]  # Valor inicial de la solución aproximada

    # Inicializar listas para almacenar valores de cada iteración
    n_vals = [c]
    xn_vals = [x0]
    fm_vals = [fe]
    dfm_vals = [dfe]
    E_vals = [error]

    # Mostrar valores iniciales
    print(f"Iteración {c}: xn = {x0}, f(xn) = {fe}, f'(xn) = {dfe}, Error = {error}")

    # Iteraciones del método de Newton
    while error > Tol and c < niter:
        xn.append(x0 - fe / dfe)  # Siguiente aproximación
        fm.append(float(f.subs(x, xn[c + 1])))  # Valor de la función en xn
        fe = fm[c + 1]  # Actualizar el valor de f
        dfm.append(float(df.subs(x, xn[c + 1])))  # Derivada en xn
        dfe = dfm[c + 1]  # Actualizar derivada

        # Cálculo del error relativo
        if abs(xn[c + 1]) > sp.S.Zero:
            E.append(abs(xn[c + 1] - x0) / abs(xn[c + 1]))
        else:
            E.append(abs(xn[c + 1] - x0))  # Fallback para evitar división por cero

        error = E[c + 1]  # Actualizar el valor del error
        x0 = xn[c + 1]  # Actualizar xn
        c += 1  # Incrementar contador

        # Almacenar valores de la iteración actual
        n_vals.append(c)
        xn_vals.append(x0)
        fm_vals.append(fe)
        dfm_vals.append(dfe)
        E_vals.append(error)

        # Mostrar valores de cada iteración
        print(f"Iteración {c}: xn = {x0}, f(xn) = {fe}, f'(xn) = {dfe}, Error = {error}")

    # Verificar condición de parada
    if fe == 0:
        print(f"{x0} es raíz de f(x)")
    elif error < Tol:
        print(f"{x0} es una aproximación de una raíz de f(x) con una tolerancia = {Tol}")
    elif dfe == 0:
        print(f"{x0} es una posible raíz múltiple de f(x)")
    else:
        print(f"Fracasó en {niter} iteraciones")

    # Crear la tabla con pandas
    T = pd.DataFrame({
        'n': n_vals,
        'xn': xn_vals,
        'fm': fm_vals,
        'E': E_vals
    })

    # Mostrar la tabla final
    print(T)
    return T
