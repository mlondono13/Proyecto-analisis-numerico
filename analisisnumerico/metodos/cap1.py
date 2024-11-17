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
    iteraciones = 0
    errores = []
    tabla = []
    while iteraciones < niter:
        fx = f.subs(x, xi)
        fpx = fprima.subs(x, xi)
        if fpx == 0:
            return None, None, None, tabla
        xi = xi - fx / fpx
        errores.append(abs(f.subs(x, xi)))
        tabla.append([iteraciones + 1, xi, abs(f.subs(x, xi))])
        iteraciones += 1
    return xi, errores, iteraciones, tabla




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
            return None, None, None, tabla
        xm = xs - fx2 * (xi - xs) / (fx1 - fx2)
        errores.append(abs(f.subs(x, xm)))
        tabla.append([iteraciones + 1, xm, abs(f.subs(x, xm))])
        if abs(f.subs(x, xm)) < Tol:
            return xm, errores, iteraciones, tabla
        xi, xs = xs, xm
        iteraciones += 1
    return xm, errores, iteraciones, tabla

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
    iteraciones_lista = list(range(1, iteraciones + 1))
    
    plt.plot(iteraciones_lista, errores, marker='o')
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
    # Lógica para procesar la interpolación o los cálculos
    resultado, errores, iteraciones, tabla = metodo_numerico(metodo, xi, xs, Tol, niter, f)
    
    df = pd.DataFrame(tabla, columns=["Iteración", "Raíz Aproximada", "Error Absoluto"])

    # Generar la gráfica
    graph_data, file_path = generar_grafica(errores, iteraciones, metodo, save_to_file=True)
    
    return resultado, df, graph_data, file_path

