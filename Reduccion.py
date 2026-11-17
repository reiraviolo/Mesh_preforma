""" Codigo que genera la matriz de puntos de una reducción """

from funciones import FuncionesGeometria
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


D1_reduccion = int(input("Ingrese el diámetro inicial del caño de la reducción: "))  # diámetro inicial del caño de la reducción - prueba(100)
D2_reduccion = int(input("Ingrese el diámetro final del caño de la reducción: "))   # diámetro final del caño de la reducción - prueba(50)
L_reduccion = int(input("Ingrese la longitud de la reducción: "))  # longitud de la reducción - prueba(150)
N_divisiones_reduccion = int(input("Ingrese el número de divisiones circunferenciales: "))  # divisiones circunferenciales - prueba(25)

nodos_reduccion, elementos_reduccion = FuncionesGeometria.Reduccion(D1_reduccion, D2_reduccion, L_reduccion, N_divisiones_reduccion)

# conectividad como enteros y en orden perimetral
elementos_reduccion_plot = elementos_reduccion.astype(int)[:, [0, 1, 3, 2]]

# coordenadas de cada vértice de cada cara: (N_elementos, 4, 3)
caras_reduccion = nodos_reduccion[elementos_reduccion_plot]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

coleccion = Poly3DCollection(
    caras_reduccion,
    facecolors="lightsteelblue",
    edgecolors="k",
    linewidths=0.6,
    alpha=0.9,
)
ax.add_collection3d(coleccion)

# nodos (opcional)
ax.scatter(nodos_reduccion[:, 0], nodos_reduccion[:, 1], nodos_reduccion[:, 2], color="r", s=10)

# límites y proporción (misma escala en los tres ejes)
mins = nodos_reduccion.min(axis=0)
maxs = nodos_reduccion.max(axis=0)
ax.set_xlim(mins[0], maxs[0])
ax.set_ylim(mins[1], maxs[1])
ax.set_zlim(mins[2], maxs[2])
ax.set_box_aspect(maxs - mins)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()