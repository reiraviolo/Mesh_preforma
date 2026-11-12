""" Codigo que genera la matriz de puntos de un codo """

from funciones import FuncionesGeometria
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


D_codo = 200  # diámetro del caño del codo
R_codo = 300  # radio de curvatura del codo (medido al eje del caño)
angulo_codo = 90  # ángulo de giro del codo [grados]
N_divisiones_codo = 25  # divisiones circunferenciales

nodos_codo, elementos_codo = FuncionesGeometria.Codo(D_codo, R_codo, angulo_codo, N_divisiones_codo)

elementos_codo_plot = elementos_codo.astype(int)[:, [0, 1, 3, 2]]

caras_codo = nodos_codo[elementos_codo_plot]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

coleccion = Poly3DCollection(
    caras_codo,
    facecolors="lightsteelblue",
    edgecolors="k",
    linewidths=0.6,
    alpha=0.9,
)
ax.add_collection3d(coleccion)

# nodos
ax.scatter(nodos_codo[:, 0], nodos_codo[:, 1], nodos_codo[:, 2], color="r", s=10)

# límites y proporción (misma escala en los tres ejes)
mins = nodos_codo.min(axis=0)
maxs = nodos_codo.max(axis=0)
ax.set_xlim(mins[0], maxs[0])
ax.set_ylim(mins[1], maxs[1])
ax.set_zlim(mins[2], maxs[2])
ax.set_box_aspect(maxs - mins)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()