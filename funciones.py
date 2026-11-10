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
