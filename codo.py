""" Codigo que genera la matriz de puntos de un codo """

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from funciones import FuncionesGeometria

if __name__ == "__main__":
    D_codo = int(input("Ingrese el diámetro del caño del codo: "))  # diámetro del caño del codo - prueba(200)
    R_codo = int(input("Ingrese el radio de curvatura del codo: "))  # radio de curvatura del codo (medido al eje del cañeo) - prueba(300)
    angulo_codo = int(input("Ingrese el ángulo de giro del codo: "))  # ángulo de giro del codo [grados] - prueba(90)
    N_divisiones_codo = int(input("Ingrese el número de divisiones circunferenciales: "))  # divisiones circunferenciales - prueba(25)

nodos_codo, elementos_codo = FuncionesGeometria.Codo(D_codo, R_codo, angulo_codo, N_divisiones_codo)

