import numpy as np

def Newtonint(x, y):
    n = len(x)
    # Inicializa la tabla con ceros, con n filas y n+1 columnas
    Tabla = np.zeros((n, n+1))
    
    # Asigna los valores de x a la primera columna y de y a la segunda columna
    Tabla[:, 0] = x
    Tabla[:, 1] = y

    # Calcula las diferencias divididas
    for j in range(2, n + 1):
        for i in range(j - 1, n):
            Tabla[i, j] = (Tabla[i, j - 1] - Tabla[i - 1, j - 1]) / (Tabla[i, 0] - Tabla[i - j + 2, 0])

    return Tabla
