import sympy as sp

def regla_falsa(f, a, b, max_iter, tol, tipErr):
    x = sp.symbols('x')
    f_sym = sp.sympify(f)
    iter = 0
    an = []
    xn = []
    bn = []
    fn = []
    E = []
    mes = ""
    err = ""

    # Verificación de entrada
    if max_iter < 0:
        err = 'El número de iteraciones es < 0'
        return iter, an, xn, bn, fn, E, mes, err
    if sp.im(f_sym.subs(x, a)):
        err = 'a no está definido en la función'
        return iter, an, xn, bn, fn, E, mes, err
    if sp.im(f_sym.subs(x, b)):
        err = 'b no está definido en la función'
        return iter, an, xn, bn, fn, E, mes, err
    if a >= b:
        err = 'a tiene que ser menor que b'
        return iter, an, xn, bn, fn, E, mes, err
    if tol < 0:
        err = 'La tolerancia es menor que 0'
        return iter, an, xn, bn, fn, E, mes, err

    # Inicialización de valores
    iter = 1
    fa = f_sym.subs(x, a).evalf()
    fb = f_sym.subs(x, b).evalf()
    m = b - (fb * (b - a)) / (fb - fa)
    fm = f_sym.subs(x, m).evalf()
    numberError = tol + 1

    # Bucle de iteración
    while numberError > tol and iter < max_iter:
        an.append(a)
        xn.append(m)
        bn.append(b)
        fn.append(fm)
        E.append(numberError)

        if sp.im(fa):
            err = 'f(a) no está definido en el dominio de la función'
            return iter, an, xn, bn, fn, E, mes, err
        if sp.im(fb):
            err = 'f(b) no está definido en el dominio de la función'
            return iter, an, xn, bn, fn, E, mes, err

        # Actualizar los valores de los puntos
        if fa * fm < 0:
            b = m
            fb = fm
        elif fm * fb < 0:
            a = m
            fa = fm

        temp = m
        if fb - fa == 0:
            err = 'La función es constante en el intervalo'
            return iter, an, xn, bn, fn, E, mes, err

        m = b - (fb * (b - a)) / (fb - fa)
        fm = f_sym.subs(x, m).evalf()
        if tipErr == 1:  # Error absoluto
            numberError = abs(m - temp)
        elif tipErr == 0:  # Error relativo
            numberError = abs((m - temp) / m)

        iter += 1

    # Almacenar la última iteración
    an.append(a)
    xn.append(m)
    bn.append(b)
    fn.append(fm)
    E.append(numberError)

    # Mensaje final
    if numberError < tol:
        mes = f"Una aproximación de la raíz fue encontrada para m = {m:.15f}"
    elif numberError > tol:
        mes = "Dado el número de iteraciones y de tolerancia, fue imposible encontrar una raíz apropiada"
    else:
        mes = "El método explotó"

    return iter, an, xn, bn, fn, E, mes, err