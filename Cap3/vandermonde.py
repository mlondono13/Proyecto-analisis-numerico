import numpy as np

def Vandermonde(x, y):
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    # Crear la matriz de Vandermonde
    A = np.vander(x)

    # Resolver el sistema lineal A * polinomio = y
    b = y.reshape(-1, 1)
    polinomio = np.linalg.solve(A, b).flatten()

    return polinomio, A
