# Proyecto-analisis-numerico

Esta es una aplicación web que permite a los usuarios interactuar con diferentes métodos numéricos para resolver problemas matemáticos. Está diseñada con una interfaz sencilla y fácil de usar. Los usuarios pueden seleccionar entre tres capítulos que cubren varios métodos numéricos, visualizar resultados y descargar gráficos generados.

## Funcionalidades

La aplicación incluye los siguientes métodos numéricos organizados en tres capítulos:

### Capítulo 1: Métodos de Bisección, Regla Falsa, Punto Fijo, Newton y Secante
- Resolución de ecuaciones no lineales mediante diferentes métodos.
- Generación de gráficos y tablas con los resultados.
- Permite ingresar funciones algebraicas y sus derivadas.

### Capítulo 2: Métodos de Jacobi, Gauss-Seidel y SOR
- Solución de sistemas de ecuaciones lineales.
- Impresión de tablas con el progreso de la solución y el radio espectral.
- Verificación de la convergencia de los métodos.

### Capítulo 3: Métodos de Interpolación (Vandermonde, Newton Interpolante, Lagrange, Spline Lineal y Cúbico)
- Interpolación de datos mediante polinomios y splines.
- Visualización de gráficos interactivos de la interpolación.
- Generación y explicación de los coeficientes de los splines (lineales y cúbicos).

## Requisitos

- **Python 3.x**: La aplicación está desarrollada con Python 3.
- **Django**: El backend está construido utilizando el framework Django.
- **Bibliotecas adicionales**:
    - `numpy`: Para cálculos numéricos.
    - `matplotlib`: Para la generación de gráficos.
    - `scipy`: Para algunos métodos avanzados.
    - `base64`: Para la conversión de gráficos a formato base64.

## Cómo Descargar y Ejecutar la Aplicación

### 1. Clonar el Repositorio

Puedes clonar el repositorio en tu máquina local utilizando Git:

```bash
git clone https://github.com/mlondono13/Proyecto-analisis-numerico.git
```

### 2. Instalar dependencias

Puedes instalarlas usando pip:

```bash
pip install requirements.txt
```

### 3. Ejecutar el servidor

Para ejecutar el servidor, ejecuta este comando en la terminal:

```bash
python manage.py runserver
```

## Uso de la Aplicación
1. En la página principal, selecciona el capítulo que deseas explorar (Capítulo 1, Capítulo 2 o Capítulo 3).
2. Ingresa los datos necesarios según el capítulo seleccionado.
3. La aplicación mostrará los resultados, gráficos y tablas correspondientes.
4. Puedes descargar los gráficos generados en formato SVG.
