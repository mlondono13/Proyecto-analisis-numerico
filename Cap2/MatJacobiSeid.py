import numpy as np

def MatJacobiSeid(x0, A, b, Tol, niter, met):
    c = 0
    error = Tol + 1
    D = np.diag(np.diag(A))
    L = -np.tril(A, -1)
    U = -np.triu(A, 1)
    E = []  # Lista para almacenar los errores en cada iteración
    
    while error > Tol and c < niter:
        if met == 0:  # Método de Jacobi
            T = np.linalg.inv(D) @ (L + U)
            C = np.linalg.inv(D) @ b
            x1 = T @ x0 + C
        elif met == 1:  # Método de Gauss-Seidel
            T = np.linalg.inv(D - L) @ U
            C = np.linalg.inv(D - L) @ b
            x1 = T @ x0 + C
        else:
            raise ValueError("El parámetro 'met' debe ser 0 (Jacobi) o 1 (Gauss-Seidel)")
        
        # Cálculo del radio espectral de T
        radio_espectral = max(abs(np.linalg.eigvals(T)))
        print(f"Radio espectral: {radio_espectral}")
        
        # Calcular el error actual
        current_error = np.linalg.norm(x1 - x0, np.inf)
        E.append(current_error)
        error = current_error
        x0 = x1
        c += 1

    # Verificación de condición de parada
    if error < Tol:
        s = x0
        print(f"Es una aproximación de la solución del sistema con una tolerancia = {Tol}")
    else:
        s = x0
        print(f"Fracasó en {niter} iteraciones")
    
    return E, s
