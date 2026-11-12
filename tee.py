import numpy as np
from funciones import FuncionesGeometria as f


OD_1 = 200  # OD straight
OD_2 = 100  # OD derivación
H_1=200 # de entrada a intersección
H_2=400 # de intersección a salida lateral
H_3=350 # de intersección a salida final
Alpha_2 = 90 # ángulo de  derivación
L_directo = H_1+H_3  # longitud de tee
L_salida = H_2  # longitud de derivación
N_directo = 25
N_salida = 25

# estas son las divisiones circunferenciaes que voy a tener en el cilindro

nodos_d, elementos_d = f.Tube(OD_1, L_directo, N_directo)
nodos_s, elementos_s = f.Tube(OD_2, L_salida, N_salida)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# conectividad como enteros y en orden perimetral
elementos = elementos.astype(int)[:, [0, 1, 3, 2]]

# coordenadas de cada vértice de cada cara: (N_elementos, 4, 3)
caras = nodos[elementos]

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

# límites y proporción
R = D_cilindro / 2
ax.set_xlim(-R, R)
ax.set_ylim(-R, R)
ax.set_zlim(0, nodos[:, 2].max())
ax.set_box_aspect((2 * R, 2 * R, nodos[:, 2].max()))

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()