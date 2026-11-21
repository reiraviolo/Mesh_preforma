"""
tubo.py — malla estructurada de un tramo de caño recto (cilindro).

Expone generar_malla_tubo(...), que a partir de las variables de entrada
devuelve la malla (nodos, elementos, grupo) y la tabla de divisiones usada.
No grafica nada: para visualizar el resultado usar graficador.py.
"""

import numpy as np


def _malla_cilindro(D, L, N, x0=0.0):
    """Malla estructurada de un cilindro de diámetro D y largo L, con eje
    en x a partir de x0. Devuelve (nodos, elementos, N_filas, paso_L) con
    elementos ya en orden perimetral [a, b, c, d]."""
    delta = np.pi * D / N
    N_segmentos_L = max(1, int(round(L / delta)))
    N_filas = N_segmentos_L + 1
    paso_L = L / N_segmentos_L

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
        x = x0 + (node_id // N) * paso_L
        y = (D / 2) * np.cos(angulo)
        z = (D / 2) * np.sin(angulo)
        nodos[i] = [x, y, z]

    return nodos, elementos, N_filas, paso_L


def generar_malla_tubo(OD, L, N=40, verbose=True):
    """
    Genera la malla de un tramo de caño recto.

    Parámetros
    ----------
    OD : float
        Diámetro exterior del caño.
    L : float
        Largo del tramo.
    N : int
        Divisiones circunferenciales (por defecto 40).

    Devuelve
    --------
    dict con:
        nodos, elementos : arrays de la malla.
        grupo : array de ceros (un solo grupo, "tubo").
        divisiones : dict con OD, L, N, N_filas (anillos), paso_L (paso axial).
    """
    nodos, elementos, N_filas, paso_L = _malla_cilindro(OD, L, N)
    grupo = np.zeros(len(elementos), dtype=int)

    divisiones = dict(OD=OD, L=L, N=N, N_filas=N_filas, paso_L=paso_L)

    if verbose:
        print(
            f"Tubo: OD={OD}, L={L}  ->  {len(nodos)} nodos, {len(elementos)} "
            f"elementos ({N} x {N_filas - 1})"
        )

    return dict(nodos=nodos, elementos=elementos, grupo=grupo, divisiones=divisiones)


if __name__ == "__main__":
    resultado = generar_malla_tubo(OD=200, L=1000, N=40)
    print(resultado["divisiones"])
