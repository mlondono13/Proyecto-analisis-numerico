import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sympy as sp
from io import BytesIO
import base64
import pandas as pd

GRAPHICS_DIR = os.path.join(os.path.dirname(__file__), 'graficas')

def generar_grafica(x, y, metodo, polinomio=None, d=None, save_to_file=False):
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, 'ro', label='Datos Originales', markersize=8)
    
    if metodo == "Spline" and polinomio is not None:
        x_vals = np.linspace(min(x), max(x), 1000)
        y_vals = np.zeros_like(x_vals)

        if d==3:
            for coef in polinomio:
                a, b, c, d, x_i = coef
                indices = (x_vals >= x_i) & (x_vals <= x_i + (x[1] - x[0]))
                x_segment = x_vals[indices] - x_i
                y_vals[indices] = a + b * x_segment + c * x_segment**2 + d * x_segment**3
        elif d==1:
            for i in range(len(x) - 1):
                x_i, x_next = x[i], x[i + 1]
                y_i, y_next = y[i], y[i + 1]

                m = (y_next - y_i) / (x_next - x_i)
                b = y_i - m * x_i
                indices = (x_vals >= x_i) & (x_vals <= x_next)
                y_vals[indices] = m * x_vals[indices] + b

        plt.plot(x_vals, y_vals, 'b-', label=f'Spline {d}-lineal' if d == 1 else 'Spline Cúbico')
    elif polinomio is not None:
        x_vals = np.linspace(min(x), max(x), 1000)
        y_vals = np.polyval(polinomio, x_vals)
        plt.plot(x_vals, y_vals, 'b-', label=f'Polinomio de Interpolación ({metodo})')
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Interpolación: {metodo}')
    plt.grid(True)
    plt.legend()

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

def Vandermonde(x, y):
    A = np.vander(x, increasing=True)
    polinomio = np.linalg.solve(A, y)
    return polinomio

def Newtonint(x, y):
    n = len(x)
    Tabla = np.zeros((n, n+1))
    Tabla[:, 0] = x
    Tabla[:, 1] = y
    for j in range(2, n+1):
        for i in range(j-1, n):
            Tabla[i, j] = (Tabla[i, j-1] - Tabla[i-1, j-1]) / (x[i] - x[i-j+1])

    coeficientes = Tabla[:, -1]
    return coeficientes

def Lagrange(x, y):
    n = len(x)
    Tabla = np.zeros((n, n))
    for i in range(n):
        Li = np.array([1])
        den = 1
        for j in range(n):
            if j != i:
                paux = np.array([1, -x[j]])
                Li = np.convolve(Li, paux)
                den *= (x[i] - x[j])
        Tabla[i, :] = y[i] * Li / den
    pol = np.sum(Tabla, axis=0)
    return pol

def Spline(x, y, d=3):
    if d == 1:
        n = len(x)
        coeficientes = []
        for i in range(n - 1):
            m = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            b = y[i] - m * x[i] 
            coeficientes.append((m, b))
        return coeficientes
    elif d == 3:
        n = len(x)
        h = [x[i + 1] - x[i] for i in range(n - 1)]
        alpha = [0] * n
        for i in range(1, n - 1):
            alpha[i] = (3 / h[i] * (y[i + 1] - y[i]) - 3 / h[i - 1] * (y[i] - y[i - 1]))
        l = [1] + [0] * (n - 1)
        mu = [0] * n
        z = [0] * n
        for i in range(1, n - 1):
            l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
            mu[i] = h[i] / l[i]
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]
        l[-1] = 1
        z[-1] = 0
        c = [0] * n
        b = [0] * (n - 1)
        d = [0] * (n - 1)
        for j in range(n - 2, -1, -1):
            c[j] = z[j] - mu[j] * c[j + 1]
            b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
            d[j] = (c[j + 1] - c[j]) / (3 * h[j])
        splines = []
        for i in range(n - 1):
            splines.append((y[i], b[i], c[i], d[i], x[i]))
        return splines
    else:
        raise ValueError("Grado no soportado. Usa d=1 (lineal) o d=3 (cúbico).")

def calcular_interpolacion(metodo, x, y, d):
    if metodo == 'Vandermonde':
        polinomio = Vandermonde(x, y)
    elif metodo == 'Newton':
        polinomio = Newtonint(x, y)
    elif metodo == 'Lagrange':
        polinomio = Lagrange(x, y)
    elif metodo == 'Spline':
        polinomio = Spline(x, y, d)
    else:
        return None

    graph_data, file_path = generar_grafica(x, y, metodo, polinomio, d, save_to_file=True)

    pol_str = f'Polinomio Resultado: {np.poly1d(polinomio)}' if metodo != 'Spline' else f'Coeficientes del Spline: {polinomio}'
    
    return pol_str, graph_data, file_path

def mostrar_resultados(metodo, x, y, d):
    pol_str, graph_data, file_path = calcular_interpolacion(metodo, x, y, d)
    return pol_str, graph_data, file_path
