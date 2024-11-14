import numpy as np

def Lagrange(x, y):
    n = len(x)
    Tabla = np.zeros((n, n))

    for i in range(n):
        Li = np.array([1])
        den = 1

        for j in range(n):
            if j != i:
                paux = np.array([1, -x[j]])
                Li = np.convolve(Li, paux)  # Convuelve los polinomios
                den *= (x[i] - x[j])  # Calcula el denominador

        Tabla[i, :] = y[i] * Li / den  # Calcula el polinomio de Lagrange para y(i)

    pol = np.sum(Tabla, axis=0)  # Suma los términos de Lagrange para obtener el polinomio final
    return pol
