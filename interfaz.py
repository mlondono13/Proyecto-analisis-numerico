import tkinter as tk
from tkinter import ttk, messagebox
from tkinter import simpledialog
from Cap1 import *
from Cap2 import *
from Cap3 import *

class NumericalMethodsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Métodos Numéricos")
        
        # Crear el menú desplegable para elegir el método
        self.method_label = ttk.Label(root, text="Seleccione el Método:")
        self.method_label.pack(pady=5)

        self.method_var = tk.StringVar()
        self.method_dropdown = ttk.Combobox(root, textvariable=self.method_var)
        self.method_dropdown['values'] = ["Lagrange", "Newton Interpolation", "Spline", "Vandermonde"]
        self.method_dropdown.pack(pady=5)
        
        # Crear botones para ingresar datos y ejecutar el método
        self.input_button = ttk.Button(root, text="Ingresar Datos", command=self.input_data)
        self.input_button.pack(pady=5)
        
        self.run_button = ttk.Button(root, text="Ejecutar Método", command=self.run_method)
        self.run_button.pack(pady=5)
        
        self.result_text = tk.Text(root, height=10, width=50)
        self.result_text.pack(pady=10)

        # Variables para almacenar los datos de entrada
        self.data_x = []
        self.data_y = []
        self.degree = 1

    def input_data(self):
        method = self.method_var.get()
        if method == "":
            messagebox.showwarning("Atención", "Seleccione un método primero.")
            return
        
        # Pedir datos específicos según el método
        if method in ["Lagrange", "Newton Interpolation", "Vandermonde"]:
            self.data_x = simpledialog.askstring("Entrada", "Ingrese valores de x separados por comas:")
            self.data_y = simpledialog.askstring("Entrada", "Ingrese valores de y separados por comas:")
            self.data_x = [float(i) for i in self.data_x.split(",")]
            self.data_y = [float(i) for i in self.data_y.split(",")]

        elif method == "Spline":
            self.data_x = simpledialog.askstring("Entrada", "Ingrese valores de x separados por comas:")
            self.data_y = simpledialog.askstring("Entrada", "Ingrese valores de y separados por comas:")
            self.degree = int(simpledialog.askstring("Entrada", "Ingrese el grado del spline (1, 2, o 3):"))
            self.data_x = [float(i) for i in self.data_x.split(",")]
            self.data_y = [float(i) for i in self.data_y.split(",")]

    def run_method(self):
        method = self.method_var.get()
        try:
            if method == "Lagrange":
                result = lagrange(self.data_x, self.data_y)
            elif method == "Newton Interpolation":
                result = newtonint(self.data_x, self.data_y)
            elif method == "Spline":
                result = spline(self.data_x, self.data_y, self.degree)
            elif method == "Vandermonde":
                result, matrix = vandermonde(self.data_x, self.data_y)
                result = f"Polinomio: {result}\nMatriz de Vandermonde:\n{matrix}"
            else:
                messagebox.showwarning("Atención", "Método no seleccionado.")
                return

            # Mostrar el resultado en el cuadro de texto
            self.result_text.delete("1.0", tk.END)
            self.result_text.insert(tk.END, f"Resultado:\n{result}")

        except Exception as e:
            messagebox.showerror("Error", f"Ocurrió un error al ejecutar el método: {e}")

if __name__ == "__main__":
    root = tk.Tk()
    app = NumericalMethodsApp(root)
    root.mainloop()