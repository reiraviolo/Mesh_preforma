import numpy as np


class FuncionesGeometria:
    @staticmethod
    def Giro(Nodos, alpha):
        # Matriz de rotación 3D alrededor del eje Y
        Rot = np.array(
            [
                [np.cos(alpha), 0, -np.sin(alpha)],
                [0, 1, 0],
                [np.sin(alpha), 0, np.cos(alpha)],
            ]
        )
        # Aplicar la rotación a cada nodo
        return np.dot(Nodos, Rot.T)

    @staticmethod
    def traslacion(Nodos, vector):
        # Aplicar la traslación a cada nodo mediante broadcasting
        return Nodos + vector

    @staticmethod
    def escalado(Nodos, factor):
        # factor puede ser un escalar o un vector (x, y, z)
        return Nodos * factor

    @staticmethod
    def Tube(D_cilindro, L_cilindro, N_divisiones):
        # Divisiones longitudinales asegurando al menos 1
        delta = np.pi * D_cilindro / N_divisiones
        N_divisiones_L = max(1, int(L_cilindro // delta))

        # MATRIZ DE NODOS (IDs)
        matriz = np.arange(N_divisiones * N_divisiones_L).reshape(
            N_divisiones_L, N_divisiones
        )

        # ARMO LA MATRIZ DE ELEMENTOS (con columna extra por ser cíclica)
        buffer = np.zeros((N_divisiones_L, N_divisiones + 1))
        buffer[:, :N_divisiones] = matriz
        buffer[:, N_divisiones] = matriz[:, 0]

        elementos = []
        for i in range(buffer.shape[0] - 1):
            for j in range(buffer.shape[1] - 1):
                elemento = buffer[i : i + 2, j : j + 2]
                elementos.append(elemento.flatten())
        elementos = np.array(elementos, dtype=int)

        # Armo la matriz de coordenadas de los nodos de forma eficiente
        buffer_nodos = matriz.flatten()
        nodos = np.zeros((len(buffer_nodos), 3))

        for i, node_id in enumerate(buffer_nodos):
            angulo = (node_id % N_divisiones) * (2 * np.pi / N_divisiones)
            x = (D_cilindro / 2) * np.cos(angulo)
            y = (D_cilindro / 2) * np.sin(angulo)
            z = (node_id // N_divisiones) * (L_cilindro / N_divisiones_L)
            nodos[i] = [x, y, z]

        return nodos, elementos

    @staticmethod
    def Codo(D_codo, R_codo, angulo_codo, N_divisiones_codo):
        theta_codo = np.deg2rad(angulo_codo) #paso de grados a radianes
        L_arco = R_codo * theta_codo
        N_divisiones_L_codo = max(1, int(L_arco // (np.pi * D_codo / N_divisiones_codo)))
        N_filas_codo = N_divisiones_L_codo + 1

        matriz_codo = np.arange(N_divisiones_codo * N_filas_codo).reshape(
            N_filas_codo, N_divisiones_codo
        )

        buffer_codo = np.zeros((N_filas_codo, N_divisiones_codo + 1))
        buffer_codo[:, :N_divisiones_codo] = matriz_codo
        buffer_codo[:, N_divisiones_codo] = matriz_codo[:, 0]

        elementos_codo = []
        for i in range(buffer_codo.shape[0] - 1):
            for j in range(buffer_codo.shape[1] - 1):
                elemento = buffer_codo[i : i + 2, j : j + 2]
                elementos_codo.append(elemento.flatten())
        elementos_codo = np.array(elementos_codo, dtype=int)

        buffer_nodos_codo = matriz_codo.flatten()
        nodos_codo = np.zeros((len(buffer_nodos_codo), 3))

        for i, node_id in enumerate(buffer_nodos_codo):
            phi = (node_id % N_divisiones_codo) * (2 * np.pi / N_divisiones_codo)
            theta = (node_id // N_divisiones_codo) * (theta_codo / N_divisiones_L_codo)
            centro = np.array([R_codo * (1 - np.cos(theta)), 0, R_codo * np.sin(theta)])
            e1 = np.array([np.cos(theta), 0, -np.sin(theta)])
            e2 = np.array([0, 1, 0])
            nodos_codo[i] = centro + (D_codo / 2) * (np.cos(phi) * e1 + np.sin(phi) * e2)

        return nodos_codo, elementos_codo