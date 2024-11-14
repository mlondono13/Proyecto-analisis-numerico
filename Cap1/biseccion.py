import sympy as sp

def biseccion(xi, xs, Tol, niter):
    x = sp.symbols('x')
    # Definir la función
    f = x**2 - 5*x + 6*sp.sin(x)
    # f = sp.exp(-x - 2) - 2*x - 2  # Alternativa de función comentada

    fi = float(f.subs(x, xi))
    fs = float(f.subs(x, xs))

    if fi == 0:
        s = xi
        E = 0
        print(f"{xi} es raíz de f(x)")
        return s, E, [fi]
    elif fs == 0:
        s = xs
        E = 0
        print(f"{xs} es raíz de f(x)")
        return s, E, [fs]
    elif fs * fi < 0:
        c = 0
        xm = (xi + xs) / 2
        fm = [float(f.subs(x, xm))]
        fe = fm[c]
        E = [Tol + 1]
        error = E[c]

        while error > Tol and fe != 0 and c < niter:
            if fi * fe < 0:
                xs = xm
                fs = float(f.subs(x, xs))
            else:
                xi = xm
                fi = float(f.subs(x, xi))

            xa = xm
            xm = (xi + xs) / 2
            fm.append(float(f.subs(x, xm)))
            fe = fm[-1]
            E.append(abs(xm - xa))
            error = E[-1]
            c += 1

        if fe == 0:
            s = xm
            print(f"{xm} es raíz de f(x)")
        elif error < Tol:
            s = xm
            print(f"{xm} es una aproximación de una raíz de f(x) con una tolerancia = {Tol}")
        else:
            s = xm
            print(f"El método falló después de {niter} iteraciones")

        return s, E, fm
    else:
        print("El intervalo es inadecuado")
        return None, None, None

# Ejemplo de uso
s, E, fm = biseccion(1, 3, 1e-5, 50)
print("Raíz aproximada:", s)
print("Errores:", E)
print("Valores de la función en cada iteración:", fm)
