import numpy as np
import pandas as pd

def MatJacobiSeid(x0, A, b, Tol, niter, met):
    c = 0
    error = Tol + 1
    D = np.diag(np.diag(A))
    L = -np.tril(A, -1)
    U = -np.triu(A, 1)
    E = []  # Lista para almacenar los errores en cada iteración
    tabla = []  # Lista para construir la tabla de resultados
    
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
        
        # Calcular el radio espectral
        radio_espectral = max(abs(np.linalg.eigvals(T)))
        
        # Guardar los resultados de la iteración
        current_error = np.linalg.norm(x1 - x0, np.inf)
        tabla.append([c + 1, *x1, current_error])
        E.append(current_error)
        
        # Actualizar error y contador
        error = current_error
        x0 = x1
        c += 1

    # Verificación de condición de parada
    convergencia = "Converge" if error < Tol else "No Converge"
    return tabla, radio_espectral, convergencia


def SOR(x0, A, b, Tol, niter, w):
    c = 0
    error = Tol + 1
    D = np.diag(np.diag(A))
    L = -np.tril(A, -1)
    U = -np.triu(A, 1)
    E = []  # Lista para almacenar los errores en cada iteración
    tabla = []  # Lista para construir la tabla de resultados
    
    while error > Tol and c < niter:
        # Calcular matrices T y C para el método SOR
        T = np.linalg.inv(D - w * L) @ ((1 - w) * D + w * U)
        C = w * np.linalg.inv(D - w * L) @ b
        x1 = T @ x0 + C
        
        # Calcular el radio espectral
        radio_espectral = max(abs(np.linalg.eigvals(T)))
        
        # Guardar los resultados de la iteración
        current_error = np.linalg.norm(x1 - x0, np.inf)
        tabla.append([c + 1, *x1, current_error])
        E.append(current_error)
        
        # Actualizar error y contador
        error = current_error
        x0 = x1
        c += 1

    # Verificación de condición de parada
    convergencia = "Converge" if error < Tol else "No Converge"
    return tabla, radio_espectral, convergencia


def mostrar_resultados_cap2(metodo, x0, A, b, Tol, niter, met=None, w=None):
    if metodo == 'jacobi' or metodo == 'gauss_seidel':
        tabla, radio_espectral, convergencia = MatJacobiSeid(x0, A, b, Tol, niter, met)
    elif metodo == 'sor':
        tabla, radio_espectral, convergencia = SOR(x0, A, b, Tol, niter, w)
    else:
        raise ValueError("Método no reconocido")
    
    # Convertir tabla a DataFrame para mejor visualización
    columns = ['Iteración'] + [f'x{i+1}' for i in range(len(x0))] + ['Error Absoluto']
    df = pd.DataFrame(tabla, columns=columns)

    print(df)
    print(f"\nRadio Espectral: {radio_espectral}")
    print(f"Resultado: {convergencia}")
    return df, radio_espectral, convergencia
