import numpy as np

def SOR(x0, A, b, Tol, niter, w):
    c = 0
    error = Tol + 1
    D = np.diag(np.diag(A))
    L = -np.tril(A, -1)
    U = -np.triu(A, 1)
    E = []  # Lista para almacenar los errores en cada iteración
    
    while error > Tol and c < niter:
        # Calcular matrices T y C para el método SOR
        T = np.linalg.inv(D - w * L) @ ((1 - w) * D + w * U)
        C = w * np.linalg.inv(D - w * L) @ b
        x1 = T @ x0 + C
        
        # Calcular el error actual
        current_error = np.linalg.norm(x1 - x0, np.inf)
        E.append(current_error)
        error = current_error
        x0 = x1
        c += 1

    # Verificación de condición de parada
    if error < Tol:
        s = x0
        p = max(abs(np.linalg.eigvals(T)))  # Radio espectral de T
        print(f"Es una aproximación de la solución del sistema con una tolerancia = {Tol}")
    else:
        s = x0
        p = None
        print(f"Fracasó en {niter} iteraciones")

    return E, s, p
