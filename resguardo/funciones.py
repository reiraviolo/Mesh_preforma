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

    @staticmethod
    def Codo(D_codo, R_codo, angulo_codo, N_divisiones_codo):
        '''Función que genera la matriz de puntos de un codo'''
        theta_codo = np.deg2rad(angulo_codo) #paso de grados a radianes
        L_arco = R_codo * theta_codo
        N_divisiones_L_codo = max(1, int(L_arco // (np.pi * D_codo / N_divisiones_codo)))
        N_filas_codo = N_divisiones_L_codo + 1

        # MATRIZ DE NODOS (IDs)
        matriz_codo = np.arange(N_divisiones_codo * N_filas_codo).reshape(
            N_filas_codo, N_divisiones_codo
        )

        # ARMO LA MATRIZ DE ELEMENTOS (con columna extra por ser cíclica)
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

    @staticmethod
    def Reduccion(D1_reduccion, D2_reduccion, L_reduccion, N_divisiones_reduccion):
        '''Función que genera la matriz de puntos de una reducción.
        Tiene 3 tramos dentro de L_reduccion:
        recto con D1 (largo 20% de D1) -> cono D1 a D2 -> recto con D2 (largo 20% de D2)'''
        L_D1 = D1_reduccion * 0.4  # un porcentaje del diámetro inicial
        L_D2 = D2_reduccion * 0.4  # un porcentaje del diámetro final

        # "Hasta" un 20%: si los tramos rectos no entran, los achico
        # proporcionalmente para que ocupen como máximo la mitad del largo
        if L_D1 + L_D2 > 0.5 * L_reduccion:
            factor = 0.5 * L_reduccion / (L_D1 + L_D2)
            L_D1 *= factor
            L_D2 *= factor
        L_cono = L_reduccion - L_D1 - L_D2

        # Divisiones longitudinales de cada tramo (al menos 1), con elementos
        # aproximadamente cuadrados según el diámetro de cada tramo
        def divisiones(largo, D):
            delta = np.pi * D / N_divisiones_reduccion
            return max(1, int(largo // delta))

        n1 = divisiones(L_D1, D1_reduccion)
        n2 = divisiones(L_cono, (D1_reduccion + D2_reduccion) / 2)
        n3 = divisiones(L_D2, D2_reduccion)

        # Posición z de cada anillo: los cambios de tramo caen justo sobre un anillo
        z1 = L_D1
        z2 = L_D1 + L_cono
        z_anillos = np.concatenate([
            np.linspace(0, z1, n1 + 1),
            np.linspace(z1, z2, n2 + 1)[1:],
            np.linspace(z2, L_reduccion, n3 + 1)[1:],
        ])
        # Diámetro de cada anillo: constante en los rectos, lineal en el cono
        D_anillos = np.interp(
            z_anillos,
            [0, z1, z2, L_reduccion],
            [D1_reduccion, D1_reduccion, D2_reduccion, D2_reduccion],
        )
        N_filas = len(z_anillos)  # = n1 + n2 + n3 + 1

        # MATRIZ DE NODOS (IDs)
        matriz = np.arange(N_divisiones_reduccion * N_filas).reshape(
            N_filas, N_divisiones_reduccion
        )

        # ARMO LA MATRIZ DE ELEMENTOS (con columna extra por ser cíclica)
        buffer = np.zeros((N_filas, N_divisiones_reduccion + 1), dtype=int)
        buffer[:, :N_divisiones_reduccion] = matriz
        buffer[:, N_divisiones_reduccion] = matriz[:, 0]

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
            angulo = (node_id % N_divisiones_reduccion) * (2 * np.pi / N_divisiones_reduccion)
            fila = node_id // N_divisiones_reduccion
            z = z_anillos[fila]
            D_actual = D_anillos[fila]
            x = (D_actual / 2) * np.cos(angulo)
            y = (D_actual / 2) * np.sin(angulo)
            nodos[i] = [x, y, z]

        return nodos, elementos