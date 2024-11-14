import numpy as np
import matplotlib.pyplot as plt

def Spline(x, y, d):
    n = len(x)
    A = np.zeros(((d + 1) * (n - 1), (d + 1) * (n - 1)))
    b = np.zeros((d + 1) * (n - 1))
    cua = x ** 2
    cub = x ** 3
    
    # Spline Lineal
    if d == 1:
        c, h = 0, 0
        for i in range(n - 1):
            A[h, c] = x[i]
            A[h, c + 1] = 1
            b[h] = y[i]
            c += 2
            h += 1
        
        c = 0
        for i in range(1, n):
            A[h, c] = x[i]
            A[h, c + 1] = 1
            b[h] = y[i]
            c += 2
            h += 1

    # Spline Cuadrático
    elif d == 2:
        c, h = 0, 0
        for i in range(n - 1):
            A[h, c] = cua[i]
            A[h, c + 1] = x[i]
            A[h, c + 2] = 1
            b[h] = y[i]
            c += 3
            h += 1
        
        c = 0
        for i in range(1, n):
            A[h, c] = cua[i]
            A[h, c + 1] = x[i]
            A[h, c + 2] = 1
            b[h] = y[i]
            c += 3
            h += 1
        
        c = 0
        for i in range(1, n - 1):
            A[h, c] = 2 * x[i]
            A[h, c + 1] = 1
            A[h, c + 3] = -2 * x[i]
            A[h, c + 4] = -1
            b[h] = 0
            c += 3
            h += 1
        
        A[h, 0] = 2
        b[h] = 0

    # Spline Cúbico
    elif d == 3:
        c, h = 0, 0
        for i in range(n - 1):
            A[h, c] = cub[i]
            A[h, c + 1] = cua[i]
            A[h, c + 2] = x[i]
            A[h, c + 3] = 1
            b[h] = y[i]
            c += 4
            h += 1

        c = 0
        for i in range(1, n):
            A[h, c] = cub[i]
            A[h, c + 1] = cua[i]
            A[h, c + 2] = x[i]
            A[h, c + 3] = 1
            b[h] = y[i]
            c += 4
            h += 1

        c = 0
        for i in range(1, n - 1):
            A[h, c] = 3 * cua[i]
            A[h, c + 1] = 2 * x[i]
            A[h, c + 2] = 1
            A[h, c + 4] = -3 * cua[i]
            A[h, c + 5] = -2 * x[i]
            A[h, c + 6] = -1
            b[h] = 0
            c += 4
            h += 1
        
        c = 0
        for i in range(1, n - 1):
            A[h, c] = 6 * x[i]
            A[h, c + 1] = 2
            A[h, c + 4] = -6 * x[i]
            A[h, c + 5] = -2
            b[h] = 0
            c += 4
            h += 1
        
        A[h, 0] = 6 * x[0]
        A[h, 1] = 2
        b[h] = 0
        h += 1
        A[h, c] = 6 * x[-1]
        A[h, c + 1] = 2
        b[h] = 0

    # Resolución del sistema
    val = np.linalg.solve(A, b)
    Tabla = val.reshape(n - 1, d + 1)

    # Gráfica de los polinomios spline
    plt.figure()
    for i in range(n - 1):
        coef = Tabla[i, :]
        x_vals = np.linspace(x[i], x[i + 1], 100)
        
        if d == 1:
            y_vals = coef[0] * x_vals + coef[1]
        elif d == 2:
            y_vals = coef[0] * x_vals**2 + coef[1] * x_vals + coef[2]
        elif d == 3:
            y_vals = coef[0] * x_vals**3 + coef[1] * x_vals**2 + coef[2] * x_vals + coef[3]
        
        plt.plot(x_vals, y_vals, 'b', linewidth=1.5)
    
    plt.plot(x, y, 'ro', markerfacecolor='r')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Spline de Interpolación (Grado {d})')
    plt.grid()
    plt.show()

    return Tabla
