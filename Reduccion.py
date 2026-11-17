""" Codigo que genera la matriz de puntos de una reducción """

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from funciones import FuncionesGeometria

if __name__ == "__main__":
    D1_reduccion = int(input("Ingrese el diámetro inicial del caño de la reducción: "))  # diámetro inicial del caño de la reducción - prueba(100)
    D2_reduccion = int(input("Ingrese el diámetro final del caño de la reducción: "))   # diámetro final del caño de la reducción - prueba(50)
    L_reduccion = int(input("Ingrese la longitud de la reducción: "))  # longitud de la reducción - prueba(150)
    N_divisiones_reduccion = int(input("Ingrese el número de divisiones circunferenciales: "))  # divisiones circunferenciales - prueba(25)

nodos_reduccion, elementos_reduccion = FuncionesGeometria.Reduccion(D1_reduccion, D2_reduccion, L_reduccion, N_divisiones_reduccion)
