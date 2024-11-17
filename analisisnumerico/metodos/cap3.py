import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sympy as sp
from io import BytesIO
import base64
import pandas as pd

# Función para generar la gráfica y convertirla a base64 para la interfaz
def generar_grafica(x, y, metodo, polinomio=None):
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, 'ro', label='Datos Originales', markersize=8)
    
    # Graficar el polinomio si se ha calculado uno
    if polinomio is not None:
        x_vals = np.linspace(min(x), max(x), 1000)
        y_vals = np.polyval(polinomio, x_vals)
        plt.plot(x_vals, y_vals, 'b-', label=f'Polinomio de Interpolación ({metodo})')
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Interpolación: {metodo}')
    plt.grid(True)
    plt.legend()

    # Guardar la gráfica como imagen SVG (vectorizable)
    img = BytesIO()
    plt.savefig(img, format='svg')
    img.seek(0)
    graph_data = base64.b64encode(img.getvalue()).decode('utf-8')  # Convertir a base64
    plt.close()  # Cerrar la figura para liberar memoria
    return graph_data

# Método de Interpolación de Vandermonde
def Vandermonde(x, y):
    A = np.vander(x, increasing=True)  # Matriz de Vandermonde
    polinomio = np.linalg.solve(A, y)  # Resolver sistema lineal A * polinomio = y
    return polinomio

# Método de Interpolación de Newton (Dividido)
def Newtonint(x, y):
    n = len(x)
    Tabla = np.zeros((n, n+1))
    Tabla[:, 0] = x
    Tabla[:, 1] = y
    for j in range(2, n+1):
        for i in range(j-1, n):
            Tabla[i, j] = (Tabla[i, j-1] - Tabla[i-1, j-1]) / (x[i] - x[i-j+1])
    
    # Los coeficientes de Newton están en la última columna de la tabla
    coeficientes = Tabla[:, -1]
    return coeficientes

# Método de Interpolación de Lagrange
def Lagrange(x, y):
    n = len(x)
    pol = np.zeros(n)
    for i in range(n):
        L = np.ones(n)
        for j in range(n):
            if j != i:
                L *= np.poly1d([1, -x[j]]) / (x[i] - x[j])  # Polinomio L(i)
        pol += y[i] * L  # Sumar al polinomio final
    return pol

# Método de Spline Lineal y Cúbico
def Spline(x, y, d=3):
    n = len(x)
    A = np.zeros((n, n))
    b = np.zeros(n)
    
    # Caso Spline Lineal (Grado 1)
    if d == 1:
        for i in range(n-1):
            A[i, i] = 1
            A[i, i+1] = x[i+1] - x[i]
            b[i] = y[i]
        A[-1, -1] = 1
        b[-1] = y[-1]
    
    # Caso Spline Cúbico (Grado 3)
    elif d == 3:
        A = np.zeros((4 * (n-1), 4 * (n-1)))
        b = np.zeros(4 * (n-1))
        for i in range(n-1):
            A[i, i] = x[i] ** 3
            A[i, i + 1] = x[i] ** 2
            A[i, i + 2] = x[i]
            A[i, i + 3] = 1
            b[i] = y[i]
        # Resolver el sistema A * coeff = b
        coeff = np.linalg.solve(A, b)
        return coeff

# Función para calcular los polinomios y generar la tabla
def calcular_interpolacion(metodo, x, y, d=1):
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
    
    # Generar la gráfica y convertirla a base64
    graph_data = generar_grafica(x, y, metodo, polinomio)
    
    # Generar el polinomio como texto para mostrar en la interfaz
    if metodo == 'Vandermonde':
        pol_str = f'Polinomio de Vandermonde: {np.poly1d(polinomio)}'
    elif metodo == 'Newton':
        pol_str = 'Polinomio de Newton: ' + ' + '.join([f'{coef}*(x-{x[i]})' for i, coef in enumerate(polinomio)])
    elif metodo == 'Lagrange':
        pol_str = 'Polinomio de Lagrange: ' + ' + '.join([f'{coef}*(x-{x[i]})' for i, coef in enumerate(polinomio)])
    else:
        pol_str = f'Coeficientes del Spline: {polinomio}'
    
    return pol_str, graph_data

# Función para mostrar los resultados en la interfaz
def mostrar_resultados(metodo, x, y, d=1):
    pol_str, graph_data = calcular_interpolacion(metodo, x, y, d)
    return pol_str, graph_data
