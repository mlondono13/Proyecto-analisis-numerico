import base64
from io import BytesIO
import os
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

GRAPHICS_DIR = os.path.join(os.path.dirname(__file__), 'graficas')

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

def punto_fijo(xi, Tol, niter, f):
    x = sp.symbols('x')
    f = sp.sympify(f)
    g = sp.solve(f - x, x)[0]
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

def newton(xi, Tol, niter, f):
    x = sp.symbols('x')
    f = sp.sympify(f)
    fprima = sp.diff(f, x)
    iteraciones = []
    errores = []
    tabla = []
    for i in range(niter):
        fx = f.subs(x, xi)
        fpx = fprima.subs(x, xi)
        if fpx == 0:
            return None, None, None, tabla

        x_next = xi - fx / fpx
        error = abs(f.subs(x, x_next))
        iteraciones.append(i + 1)
        errores.append(error)
        tabla.append([i + 1, x_next, error])

        if error < Tol:
            break
        xi = x_next
    return xi, errores, iteraciones, tabla


def secante(xi, xs, Tol, niter, f):
    x = sp.symbols('x')
    f = sp.sympify(f)
    iteraciones = []
    errores = []
    tabla = []

    for i in range(niter):
        fxi = f.subs(x, xi)
        fxs = f.subs(x, xs)
        if (fxs - fxi) == 0:
            return None, None, None, tabla
        x_next = xs - (fxs * (xs - xi)) / (fxs - fxi)
        error = abs(x_next - xs)
        iteraciones.append(i + 1)
        errores.append(error)
        tabla.append([i + 1, x_next, error])
        if error < Tol:
            break
        xi, xs = xs, x_next
    return x_next, errores, iteraciones, tabla

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

def generar_grafica(errores, iteraciones, metodo, save_to_file=False):
    plt.plot(iteraciones, errores, marker='o')
    plt.yscale('log')
    plt.xlabel('Iteración')
    plt.ylabel('Error Absoluto')
    plt.title(f'Convergencia del Método de {metodo.capitalize()}')
    plt.grid(True)

    img = BytesIO()
    plt.savefig(img, format='svg')
    
    if save_to_file:
        file_path = os.path.join(GRAPHICS_DIR, f'grafica_{metodo}.svg')
        with open(file_path, 'wb') as f:
            f.write(img.getvalue())
        img.seek(0)
        graph_data = base64.b64encode(img.getvalue()).decode('utf-8')
        return graph_data, file_path 
    else:
        img.seek(0)
        graph_data = base64.b64encode(img.getvalue()).decode('utf-8')
        plt.close()
        return graph_data, None

def mostrar_resultados1(metodo, xi, xs, Tol, niter, f):
    resultado, errores, iteraciones, tabla = metodo_numerico(metodo, xi, xs, Tol, niter, f)
    df = pd.DataFrame(tabla, columns=["Iteración", "Raíz Aproximada", "Error Absoluto"])
    graph_data, file_path = generar_grafica(errores, iteraciones, metodo, save_to_file=True)
    return resultado, df, graph_data, file_path

