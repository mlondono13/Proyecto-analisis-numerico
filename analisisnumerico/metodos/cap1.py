import base64
from io import BytesIO
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# Método de Bisección
def biseccion(xi, xs, Tol, niter, f):
    x = sp.symbols('x')
    f = sp.sympify(f)
    iteraciones = 0
    errores = []
    tabla = []
    while iteraciones < niter:
        xm = (xi + xs) / 2
        fxm = f.subs(x, xm)
        if abs(fxm) < Tol:
            tabla.append([iteraciones + 1, xm, abs(f.subs(x, xm))])
            return xm, errores, iteraciones, tabla
        if f.subs(x, xi) * fxm < 0:
            xs = xm
        else:
            xi = xm
        errores.append(abs(f.subs(x, xm)))
        tabla.append([iteraciones + 1, xm, abs(f.subs(x, xm))])
        iteraciones += 1
    return xm, errores, iteraciones, tabla

# Método de la Regla Falsa
def regla_falsa(xi, xs, Tol, niter, f):
    x = sp.symbols('x')
    f = sp.sympify(f)
    iteraciones = 0
    errores = []
    tabla = []
    while iteraciones < niter:
        xm = xs - f.subs(x, xs) * (xi - xs) / (f.subs(x, xi) - f.subs(x, xs))
        fxm = f.subs(x, xm)
        if abs(fxm) < Tol:
            tabla.append([iteraciones + 1, xm, abs(f.subs(x, xm))])
            return xm, errores, iteraciones, tabla
        if f.subs(x, xi) * fxm < 0:
            xs = xm
        else:
            xi = xm
        errores.append(abs(f.subs(x, xm)))
        tabla.append([iteraciones + 1, xm, abs(f.subs(x, xm))])
        iteraciones += 1
    return xm, errores, iteraciones, tabla

# Método de Punto Fijo
def punto_fijo(xi, Tol, niter, f):
    x = sp.symbols('x')
    f = sp.sympify(f)
    g = sp.solve(f - x, x)[0]  # Solucionamos f(x) = x para obtener g(x)
    iteraciones = 0
    errores = []
    tabla = []
    while iteraciones < niter:
        xi = g.subs(x, xi)
        errores.append(abs(f.subs(x, xi)))
        tabla.append([iteraciones + 1, xi, abs(f.subs(x, xi))])
        if abs(f.subs(x, xi)) < Tol:
            return xi, errores, iteraciones, tabla
        iteraciones += 1
    return xi, errores, iteraciones, tabla

# Método de Newton-Raphson
def newton(xi, Tol, niter, f):
    x = sp.symbols('x')
    f = sp.sympify(f)
    fprima = sp.diff(f, x)  # Derivada de f
    iteraciones = 0
    errores = []
    tabla = []
    while iteraciones < niter:
        fx = f.subs(x, xi)
        fpx = fprima.subs(x, xi)
        if fpx == 0:
            return None, None, None, tabla  # Evitar división por cero
        xi = xi - fx / fpx
        errores.append(abs(f.subs(x, xi)))
        tabla.append([iteraciones + 1, xi, abs(f.subs(x, xi))])
        if abs(f.subs(x, xi)) < Tol:
            return xi, errores, iteraciones, tabla
        iteraciones += 1
    return xi, errores, iteraciones, tabla

# Método de la Secante
def secante(xi, xs, Tol, niter, f):
    x = sp.symbols('x')
    f = sp.sympify(f)
    iteraciones = 0
    errores = []
    tabla = []
    while iteraciones < niter:
        fx1 = f.subs(x, xi)
        fx2 = f.subs(x, xs)
        if fx2 - fx1 == 0:
            return None, None, None, tabla  # Evitar división por cero
        xm = xs - fx2 * (xi - xs) / (fx1 - fx2)
        errores.append(abs(f.subs(x, xm)))
        tabla.append([iteraciones + 1, xm, abs(f.subs(x, xm))])
        if abs(f.subs(x, xm)) < Tol:
            return xm, errores, iteraciones, tabla
        xi, xs = xs, xm
        iteraciones += 1
    return xm, errores, iteraciones, tabla

# Función controladora
def metodo_numerico(metodo, xi, xs, Tol, niter, f):
    x = sp.symbols('x')
    f = sp.sympify(f)

    if metodo == 'biseccion':
        return biseccion(xi, xs, Tol, niter, f)
    elif metodo == 'regla_falsa':
        return regla_falsa(xi, xs, Tol, niter, f)
    elif metodo == 'punto_fijo':
        return punto_fijo(xi, Tol, niter, f)
    elif metodo == 'newton':
        return newton(xi, Tol, niter, f)
    elif metodo == 'secante':
        return secante(xi, xs, Tol, niter, f)
    else:
        return None, None, None, None

def generar_grafica(errores, iteraciones, metodo):
    plt.plot(range(1, iteraciones + 1), errores, marker='o')
    plt.yscale('log')  # Escala logarítmica
    plt.xlabel('Iteración')
    plt.ylabel('Error Absoluto')
    plt.title(f'Convergencia del Método de {metodo.capitalize()}')
    plt.grid(True)

    # Guardar la imagen en formato SVG en lugar de PNG
    img = BytesIO()
    plt.savefig(img, format='svg')
    img.seek(0)
    graph_data = base64.b64encode(img.getvalue()).decode('utf-8')  # Convertir a base64
    plt.close()  # Cerrar la figura para liberar memoria
    return graph_data

# Mostrar la tabla y gráfica
def mostrar_resultados(metodo, xi, xs, Tol, niter, f):
    resultado, errores, iteraciones, tabla = metodo_numerico(metodo, xi, xs, Tol, niter, f)
    
    # Crear DataFrame para la tabla
    df = pd.DataFrame(tabla, columns=["Iteración", "Raíz Aproximada", "Error Absoluto"])
    print(df)
    
    # Generar la gráfica en base64 en formato SVG
    graph_data = generar_grafica(errores, iteraciones, metodo)
    
    return resultado, errores, iteraciones, tabla, df, graph_data
