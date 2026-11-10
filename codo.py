""" Codigo que genera la matriz de puntos de un codo """
import numpy as np
import matplotlib as plt


D_codo = 200  # diámetro del caño del codo
R_codo = 300  # radio de curvatura del codo (medido al eje del caño)
angulo_codo = 90  # ángulo de giro del codo [grados]
N_divisiones_codo = 25  # divisiones circunferenciales

theta_codo = np.deg2rad(angulo_codo)
# divisiones a lo largo del codo: igual que en el cilindro, busco elementos
# aproximadamente cuadrados usando la longitud del arco sobre el eje
L_arco = R_codo * theta_codo
N_divisiones_L_codo = int(L_arco // (np.pi * D_codo / N_divisiones_codo))
# filas de nodos = divisiones + 1 (para que el codo cubra el ángulo completo)
N_filas_codo = N_divisiones_L_codo + 1

# MATRIZ DE NODOS
matriz_codo = np.arange(N_divisiones_codo * N_filas_codo).reshape(
    N_filas_codo, N_divisiones_codo
)

# ARMO LA MATRIZ DE ELEMENTOS
# creo la matriz de ceros + una columna por ser cíclica
buffer_codo = np.zeros((N_filas_codo, N_divisiones_codo + 1))
# ingreso los valores de la matriz original
buffer_codo[:, 0:N_divisiones_codo] = matriz_codo
# ingreso los valores a la última columna del buffer
buffer_codo[:, N_divisiones_codo] = matriz_codo[:, 0]

elementos_codo = []
for i in range(buffer_codo.shape[0] - 1):
    for j in range(buffer_codo.shape[1] - 1):
        elemento = buffer_codo[i : i + 2, j : j + 2]
        elementos_codo.append(elemento.flatten())
elementos_codo = np.array(elementos_codo)

# Armo la matriz de coordenadas de los nodos
# El eje del codo arranca en el origen apuntando en +z y gira hacia +x,
# con centro de curvatura en (R_codo, 0, 0). Así la sección en theta = 0
# coincide con la base del cilindro (plano z = 0).
r = D_codo / 2
buffer_nodos_codo = matriz_codo.flatten()
nodos_codo = np.zeros((len(buffer_nodos_codo), 3))
for i in range(len(buffer_nodos_codo)):
    # ángulo alrededor del caño
    phi = (buffer_nodos_codo[i] % N_divisiones_codo) * (2 * np.pi / N_divisiones_codo)
    # ángulo a lo largo del codo
    theta = (buffer_nodos_codo[i] // N_divisiones_codo) * (theta_codo / N_divisiones_L_codo)
    # punto sobre el eje del codo
    centro = np.array([R_codo * (1 - np.cos(theta)), 0, R_codo * np.sin(theta)])
    # vectores de la sección transversal (perpendiculares a la tangente)
    e1 = np.array([np.cos(theta), 0, -np.sin(theta)])
    e2 = np.array([0, 1, 0])
    nodos_codo[i] = centro + r * np.cos(phi) * e1 + r * np.sin(phi) * e2

#print(nodos_codo)