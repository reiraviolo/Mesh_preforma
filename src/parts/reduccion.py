"""
reduccion.py — malla estructurada de una reducción (cono) entre dos diámetros.

Expone generar_malla_reduccion(...), que a partir de las variables de entrada
devuelve la malla (nodos, elementos, grupo) y la tabla de divisiones usada.
No grafica nada: para visualizar el resultado usar graficador.py.
"""

import numpy as np


def _malla_reduccion(OD_1, OD_2, L, N):
    """Malla estructurada de una reducción, con eje en X (mismo eje de
    entrada que tubo.py/codo.py/tee.py), entre OD_1 (x=0) y OD_2 (x=L).
    Tiene 3 tramos: recto con OD_1 -> cono OD_1 a OD_2 -> recto con OD_2.
    Devuelve (nodos, elementos, info) con elementos ya en orden perimetral
    [a, b, c, d]."""
    L_D1 = OD_1 * 0.4
    L_D2 = OD_2 * 0.4

    # Si los tramos rectos no entran, los achico para que ocupen como
    # máximo la mitad del largo total.
    if L_D1 + L_D2 > 0.5 * L:
        factor = 0.5 * L / (L_D1 + L_D2)
        L_D1 *= factor
        L_D2 *= factor
    L_cono = L - L_D1 - L_D2

    def divisiones_tramo(largo, D):
        delta = np.pi * D / N
        return max(1, int(largo // delta))

    n1 = divisiones_tramo(L_D1, OD_1)
    n2 = divisiones_tramo(L_cono, (OD_1 + OD_2) / 2)
    n3 = divisiones_tramo(L_D2, OD_2)

    x1 = L_D1
    x2 = L_D1 + L_cono
    x_anillos = np.concatenate(
        [
            np.linspace(0, x1, n1 + 1),
            np.linspace(x1, x2, n2 + 1)[1:],
            np.linspace(x2, L, n3 + 1)[1:],
        ]
    )
    D_anillos = np.interp(x_anillos, [0, x1, x2, L], [OD_1, OD_1, OD_2, OD_2])
    N_filas = len(x_anillos)

    matriz = np.arange(N * N_filas).reshape(N_filas, N)
    buffer = np.zeros((N_filas, N + 1), dtype=int)
    buffer[:, :N] = matriz
    buffer[:, N] = matriz[:, 0]

    elementos = []
    for i in range(buffer.shape[0] - 1):
        for j in range(buffer.shape[1] - 1):
            elementos.append(buffer[i : i + 2, j : j + 2].flatten())
    elementos = np.array(elementos, dtype=int)[:, [0, 1, 3, 2]]

    ids = matriz.flatten()
    nodos = np.zeros((len(ids), 3))
    for i, node_id in enumerate(ids):
        angulo = (node_id % N) * (2 * np.pi / N)
        fila = node_id // N
        D_actual = D_anillos[fila]
        x = x_anillos[fila]
        y = (D_actual / 2) * np.cos(angulo)
        z = (D_actual / 2) * np.sin(angulo)
        nodos[i] = [x, y, z]

    info = dict(n1=n1, n2=n2, n3=n3, N_filas=N_filas, L_D1=L_D1, L_cono=L_cono, L_D2=L_D2)
    return nodos, elementos, info


def generar_malla_reduccion(OD_1, OD_2, L, N=40, verbose=True):
    """
    Genera la malla de una reducción (cono) entre dos diámetros.

    Parámetros
    ----------
    OD_1 : float
        Diámetro exterior en el extremo de entrada (x=0).
    OD_2 : float
        Diámetro exterior en el extremo de salida (x=L).
    L : float
        Largo total de la reducción.
    N : int
        Divisiones circunferenciales (por defecto 40).

    Devuelve
    --------
    dict con:
        nodos, elementos : arrays de la malla.
        grupo : array de ceros (un solo grupo, "reducción").
        divisiones : dict con OD_1, OD_2, L, N, N_filas, y el detalle de los
        3 tramos (n1/n2/n3 divisiones, L_D1/L_cono/L_D2 largos).
    """
    nodos, elementos, info = _malla_reduccion(OD_1, OD_2, L, N)
    grupo = np.zeros(len(elementos), dtype=int)

    divisiones = dict(OD_1=OD_1, OD_2=OD_2, L=L, N=N, **info)

    if verbose:
        print(
            f"Reducción: OD_1={OD_1} -> OD_2={OD_2}, L={L}  ->  {len(nodos)} "
            f"nodos, {len(elementos)} elementos"
        )

    return dict(nodos=nodos, elementos=elementos, grupo=grupo, divisiones=divisiones)


if __name__ == "__main__":
    resultado = generar_malla_reduccion(OD_1=100, OD_2=50, L=150, N=25)
    print(resultado["divisiones"])
