import sympy as sp

def raices_multiples(f, x0, tol, max_iter, tipErr):
    """
    Método de raíces múltiples para encontrar una raíz de la función `f`.
    
    Parámetros:
    f: str - Función en formato de string (ej. "pi**(-x)*(-1+x)+x**(2/3)-19").
    x0: float - Valor inicial.
    tol: float - Tolerancia.
    max_iter: int - Número máximo de iteraciones.
    tipErr: int - Tipo de error (1 para error absoluto, otro valor para error relativo).
    
    Retorna:
    c: int - Número de iteraciones realizadas.
    xn: list - Lista de aproximaciones de la raíz.
    fm: list - Valores de la función en las iteraciones.
    dfm: list - Derivada primera de la función en las iteraciones.
    d2fm: list - Derivada segunda de la función en las iteraciones.
    E: list - Lista de errores.
    mes: str - Mensaje final del resultado.
    err: str - Mensaje de error (si aplica).
    """
    # Definir símbolos y función simbólica
    x = sp.Symbol('x')
    f_sym = sp.sympify(f)
    df = sp.diff(f_sym, x)
    d2f = sp.diff(df, x)

    # Inicialización de variables
    c = 0
    xn = [x0]
    fm = [float(f_sym.subs(x, x0))]
    dfm = [float(df.subs(x, x0))]
    d2fm = [float(d2f.subs(x, x0))]
    E = [tol + 1]
    mes = ""
    err = ""

    # Validaciones iniciales
    if sp.im(fm[0]) != 0:
        err = "f(x0) no está definido en el dominio de la función"
        return c, xn, fm, dfm, d2fm, E, mes, err

    if max_iter < 0:
        err = "El número de iteraciones es menor que 0"
        return c, xn, fm, dfm, d2fm, E, mes, err

    if tol < 0:
        err = "La tolerancia es menor que 0"
        return c, xn, fm, dfm, d2fm, E, mes, err

    # Iteración del método
    while E[c] > tol and c < max_iter:
        fe = fm[c]
        dfe = dfm[c]
        d2fe = d2fm[c]
        d = dfe**2 - fe * d2fe

        if d == 0:
            err = f"División por cero en la iteración: {c}"
            return c, xn, fm, dfm, d2fm, E, mes, err

        # Calcular la siguiente aproximación
        x_next = x0 - (fe * dfe) / d
        xn.append(float(x_next))

        if x_next == float('inf') or x_next == float('-inf'):
            err = f"Valor infinito en la iteración: {c}"
            return c, xn, fm, dfm, d2fm, E, mes, err

        # Actualizar valores
        fm.append(float(f_sym.subs(x, x_next)))
        dfm.append(float(df.subs(x, x_next)))
        d2fm.append(float(d2f.subs(x, x_next)))

        if tipErr == 1:  # Error absoluto
            E.append(abs(x_next - x0))
        else:  # Error relativo
            E.append(abs(x_next - x0) / abs(x_next))

        if sp.im(fm[c + 1]) != 0 or sp.im(dfm[c + 1]) != 0 or sp.im(d2fm[c + 1]) != 0:
            err = f"xi no está definido en el dominio de la función en la iteración: {c + 1}"
            return c, xn, fm, dfm, d2fm, E, mes, err

        # Preparar para la siguiente iteración
        x0 = x_next
        c += 1

    # Mensaje final
    if fm[-1] == 0:
        mes = f"La raíz fue encontrada en x = {xn[-1]:.15f}"
    elif E[-1] <= tol:
        mes = f"Una aproximación de la raíz fue encontrada en x = {xn[-1]:.15f}"
    elif c == max_iter:
        mes = "Dado el número de iteraciones y la tolerancia, fue imposible encontrar una raíz apropiada"
    else:
        mes = "El método explotó"

    return c, xn, fm, dfm, d2fm, E, mes, err
