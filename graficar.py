"""Esta función permite graficar lo que el main requiera"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def graficar(nodos, elementos, titulo):
    """Función para graficar los nodos y elementos de un accesorio en 3D."""

    # conectividad como enteros y en orden perimetral
    elementos_plot = elementos.astype(int)[:, [0, 1, 3, 2]]

    # coordenadas de cada vértice de cada cara: (N_elementos, 4, 3)
    caras = nodos[elementos_plot]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

    coleccion = Poly3DCollection(
        caras,
        facecolors="lightsteelblue",
        edgecolors="k",
        linewidths=0.6,
        alpha=0.9,
    )
    ax.add_collection3d(coleccion)

    # nodos (opcional)
    ax.scatter(nodos[:, 0], nodos[:, 1], nodos[:, 2], color="r", s=10)

    # límites y proporción (misma escala en los tres ejes)
    mins = nodos.min(axis=0)
    maxs = nodos.max(axis=0)
    ax.set_xlim(mins[0], maxs[0])
    ax.set_ylim(mins[1], maxs[1])
    ax.set_zlim(mins[2], maxs[2])
    ax.set_box_aspect(maxs - mins)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.title(titulo)
    