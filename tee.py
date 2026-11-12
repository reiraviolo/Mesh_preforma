import numpy as np
from funciones import FuncionesGeometria as f


OD_1 = 200  # OD straight
OD_2 = 150  # OD derivación
H_1 = 200  # de entrada a intersección
H_2 = 400  # de intersección a salida lateral
H_3 = 350  # de intersección a salida final
Alpha_2 = 45  # ángulo de  derivación
L_directo = H_1 + H_3  # longitud de tee
L_salida = H_2  # longitud de derivación
N_directo = 25
N_salida = 25

# estas son las divisiones circunferenciaes que voy a tener en el cilindro

# genero tubo recto
nodos_d, elementos_d = f.Tube(OD_1, L_directo, N_directo)
nodos_d = f.traslacion(nodos_d, np.array([-H_1, 0, 0]))
# genero tubo lateral y giro para que quede en la posición correcta
nodos_s, elementos_s = f.Tube(OD_2, L_salida, N_salida)
nodos_s = f.Giro(nodos_s, np.radians(Alpha_2))
# print("paso 1")


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# conectividad como enteros y en orden perimetral
elementos_d = elementos_d.astype(int)[:, [0, 1, 3, 2]]
elementos_s = elementos_s.astype(int)[:, [0, 1, 3, 2]]

# coordenadas de cada vértice de cada cara: (N_elementos, 4, 3)
caras_d = nodos_d[elementos_d]
caras_s = nodos_s[elementos_s]

# print("paso 2")

# curva intersección
# limite de evaliación -OD_2 -> +OD_2, Condición: OD_2<=OD_1
EV = range(-OD_2 // 2, (OD_2 // 2) + 1, 1)
print(EV)
x = []
y = []
z = []

# El eje x es el eje del tubo directo, por lo que la coordenada transversal
# (la que se recorre de -r2 a +r2) ahora es y, y x es la que se calcula.
for y_val in EV:
    Rp = np.sqrt((OD_2 / 2) ** 2 - y_val**2) / np.cos(np.radians(Alpha_2))
    Z_1 = np.sqrt((OD_1 / 2) ** 2 - y_val**2)
    x.append((Z_1 - Rp) / np.tan(np.radians(Alpha_2)))
    y.append(y_val)
    z.append(Z_1)

EV = range((OD_2 // 2), (-OD_2 // 2) - 1, -1)
for y_val in EV:
    Rp = np.sqrt((OD_2 / 2) ** 2 - y_val**2) / np.cos(np.radians(Alpha_2))
    Z_1 = np.sqrt((OD_1 / 2) ** 2 - y_val**2)
    x.append((Z_1 + Rp) / np.tan(np.radians(Alpha_2)))
    y.append(y_val)
    z.append(Z_1)
#fin curva intersección

#filtro de nodos de salida para que no se vean los que están dentro del tubo directo
mask_nodos_s = nodos_s[:, 1]**2+ nodos_s[:, 2]**2 > OD_1**2/4

#filtro de nodos del directo para que no se vean los que están dentro del tubo directo
nodos_dm = f.Giro(nodos_d, np.radians(-Alpha_2))   # solo para evaluar la condición

dentro = (
    (nodos_dm[:, 1]**2 + nodos_dm[:, 2]**2 < OD_2**2/4)   # dentro del cilindro lateral
    & (nodos_dm[:, 0] > 0)                                # lado de la derivación
    & (nodos_d[:, 0] > 0)                                 # mitad positiva del tubo recto
)
mask_nodos_d = ~dentro                                    # True = el nodo se queda

# salida de la intersección
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

coleccion_s = Poly3DCollection(
    caras_s,
    facecolors="lightsteelblue",
    edgecolors="k",
    linewidths=0.6,
    alpha=0.9,
)
#ax.add_collection3d(coleccion_s)

coleccion_d = Poly3DCollection(
    caras_d,
    facecolors="lightsteelblue",
    edgecolors="k",
    linewidths=0.6,
    alpha=0.9,
)
#ax.add_collection3d(coleccion_d)

# nodos (opcional)
#ax.scatter(nodos_d[:, 0], nodos_d[:, 1], nodos_d[:, 2], color="r", s=10)
#ax.scatter(nodos_s[:, 0], nodos_s[:, 1], nodos_s[:, 2], color="r", s=10)
ax.scatter(nodos_s[mask_nodos_s, 0], nodos_s[mask_nodos_s, 1], nodos_s[mask_nodos_s, 2], color="r", s=10)
ax.scatter(nodos_d[mask_nodos_d, 0], nodos_d[mask_nodos_d, 1], nodos_d[mask_nodos_d, 2], color="r", s=10)

# Graficar la curva en 3D
ax.scatter(x, y, z, label="Curva de Intersección", color="b", linewidth=2)



# límites y proporción
todos = np.vstack([nodos_d, nodos_s])
mins = todos.min(axis=0)
maxs = todos.max(axis=0)
ax.set_xlim(mins[0], maxs[0])
ax.set_ylim(mins[1], maxs[1])
ax.set_zlim(mins[2], maxs[2])
ax.set_box_aspect(maxs - mins)  # proporción real (x = eje del tubo directo)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()
