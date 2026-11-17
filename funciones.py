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
        # Segmentos longitudinales (al menos 1), con elementos lo más cuadrados posible.
        # Hay un anillo más que segmentos, así el último anillo queda en x = L.
        delta = np.pi * D_cilindro / N_divisiones
        N_segmentos_L = max(1, int(round(L_cilindro / delta)))
        N_divisiones_L = N_segmentos_L + 1  # cantidad de anillos de nodos
        paso_L = L_cilindro / N_segmentos_L

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
            x = (node_id // N_divisiones) * paso_L
            y = (D_cilindro / 2) * np.cos(angulo)
            z = (D_cilindro / 2) * np.sin(angulo)

            nodos[i] = [x, y, z]

        return nodos, elementos
