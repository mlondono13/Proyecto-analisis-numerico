from django.shortcuts import render, redirect
import numpy as np
import pandas as pd

from .cap1 import *
from .cap2 import *
from .cap3 import *

def inicio(request):
    if request.method == "POST":
        capitulo = request.POST.get('capitulo')
        if capitulo == '1':
            return redirect('capitulo_1')
        elif capitulo == '2':
            return redirect('capitulo_2')
        elif capitulo == '3':
            return redirect('capitulo_3')
    return render(request, 'metodos/inicio.html')

def calcular_cap1(request):
    if request.method == "POST":
        metodo = request.POST.get('metodo')
        xi = float(request.POST.get('xi'))
        xs = float(request.POST.get('xs'))
        Tol = float(request.POST.get('tol'))
        niter = int(request.POST.get('niter'))
        f = request.POST.get('funcion')

        # Llamar a la función para calcular el resultado y generar la gráfica
        resultado, errores, iteraciones, tabla, df, graph_data = mostrar_resultados(metodo, xi, xs, Tol, niter, f)

        # Devolver la plantilla con los resultados y la gráfica
        return render(request, 'metodos/cap1.html', {
            'resultado': resultado,
            'df': df.to_html(classes="table table-bordered table-striped"),
            'graph_data': graph_data,
            'metodo': metodo
        })
    
    return render(request, 'metodos/cap1.html')

def calcular_cap2(request):
    resultado = None
    if request.method == 'POST':
        try:
            metodo = request.POST.get('metodo')

            # Obtener y procesar las matrices
            A = np.array([list(map(float, row.split())) for row in request.POST.get('matriz_A').split(';')])
            b = np.array(list(map(float, request.POST.get('vector_b').split(','))))
            x0 = np.array(list(map(float, request.POST.get('x0').split(','))))
            Tol = float(request.POST.get('tolerancia'))
            niter = int(request.POST.get('iteraciones'))
            w = float(request.POST.get('w')) if request.POST.get('w') else None

            # Ejecutar el método seleccionado
            if metodo == 'jacobi':
                tabla, radio_espectral, convergencia = MatJacobiSeid(x0, A, b, Tol, niter, met=0)
            elif metodo == 'gauss_seidel':
                tabla, radio_espectral, convergencia = MatJacobiSeid(x0, A, b, Tol, niter, met=1)
            elif metodo == 'sor':
                if w is None:
                    w = 1.0
                tabla, radio_espectral, convergencia = SOR(x0, A, b, Tol, niter, w)
            else:
                raise ValueError("Método no válido")

            # Convertir la tabla a un DataFrame para su visualización en HTML
            columns = ['Iteración'] + [f'x{i+1}' for i in range(len(x0))] + ['Error Absoluto']
            df = pd.DataFrame(tabla, columns=columns)

            resultado = {
                'tabla': df.to_html(classes="table table-bordered table-striped"),
                'radio_espectral': radio_espectral,
                'convergencia': convergencia,
                'metodo': metodo
            }

        except Exception as e:
            resultado = {'error': f"Error al procesar los datos: {e}"}

    return render(request, 'metodos/cap2.html', {'resultado': resultado})


def calcular_cap3(request):
    if request.method == "POST":
        metodo = request.POST.get('metodo')
        x_values = list(map(float, request.POST.get('x_values').split(',')))
        y_values = list(map(float, request.POST.get('y_values').split(',')))

        # Calcular los polinomios y generar la gráfica
        pol_str, graph_data = mostrar_resultados(metodo, x_values, y_values)

        # Devolver la plantilla con los resultados y la gráfica
        return render(request, 'metodos/cap3.html', {
            'pol_str': pol_str,
            'graph_data': graph_data,
            'metodo': metodo
        })
    
    return render(request, 'metodos/cap3.html')