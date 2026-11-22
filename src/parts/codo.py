"""
codo.py — malla estructurada de un codo (curva toroidal) de diámetro constante.

Expone generar_malla_codo(...), que a partir de las variables de entrada
devuelve la malla (nodos, elementos, grupo) y la tabla de divisiones usada.
No grafica nada: para visualizar el resultado usar graficador.py.
"""

import numpy as np


def _malla_codo(D, R, angulo, N):
    """Malla estructurada de un codo de diámetro D, radio de curvatura R
    (medido al eje del caño) y ángulo de giro `angulo` (grados). El codo
    arranca en el origen con su eje tangente a X (mismo eje de entrada que
    tubo.py/reduccion.py/tee.py), y gira en el plano X-Z alrededor del eje
    Y. Devuelve (nodos, elementos, N_filas, N_filas_L) con elementos ya en
    orden perimetral [a, b, c, d]."""
    theta = np.deg2rad(angulo)
    L_arco = R * theta
    N_filas_L = max(1, int(L_arco // (np.pi * D / N)))
    N_filas = N_filas_L + 1

    matriz = np.arange(N * N_filas).reshape(N_filas, N)
    buffer = np.zeros((N_filas, N + 1), dtype=int)
    buffer[:, :N] = matriz
    buffer[:, N] = matriz[:, 0]

    elementos = []
    for i in range(buffer.shape[0] - 1):
        for j in range(buffer.shape[1] - 1):
            elementos.append(buffer[i : i + 2, j : j + 2].flatten())
    # [2,3,1,0] en vez de [0,1,3,2]: con la nueva parametrizacion (entrada
    # en +X) ese orden recorre el perimetro al reves que tubo/reduccion/tee;
    # se verifico con control_malla (orientacion_inconsistente pasaba de
    # 48 a 0 en la red ensamblada al invertir el orden).
    elementos = np.array(elementos, dtype=int)[:, [2, 3, 1, 0]]

    ids = matriz.flatten()
    nodos = np.zeros((len(ids), 3))
    for i, node_id in enumerate(ids):
        phi = (node_id % N) * (2 * np.pi / N)
        t = (node_id // N) * (theta / N_filas_L)
        # centro(t=0) = origen, con tangente (R,0,0) -> entra por +X
        centro = np.array([R * np.sin(t), 0, R * (1 - np.cos(t))])
        # e1/e2: base del circulo de la seccion transversal, perpendicular
        # a la tangente en cada t (e1 gira junto con la curva en el plano
        # X-Z, e2 queda fijo en Y, perpendicular al plano de giro)
        e1 = np.array([-np.sin(t), 0, np.cos(t)])
        e2 = np.array([0, 1, 0])
        nodos[i] = centro + (D / 2) * (np.cos(phi) * e1 + np.sin(phi) * e2)

    return nodos, elementos, N_filas, N_filas_L


def generar_malla_codo(OD, R, angulo=90.0, N=40, verbose=True):
    """
    Genera la malla de un codo (curva de 90° u otro ángulo).

    Parámetros
    ----------
    OD : float
        Diámetro exterior del caño.
    R : float
        Radio de curvatura, medido al eje del caño.
    angulo : float
        Ángulo de giro del codo, en grados (por defecto 90).
    N : int
        Divisiones circunferenciales (por defecto 40).

    El codo arranca en el origen con tangente de entrada en +X (el mismo eje
    de entrada que tubo.py/reduccion.py/tee.py) y gira hacia +Z en el plano
    x-z, alrededor del eje Y.

    Devuelve
    --------
    dict con:
        nodos, elementos : arrays de la malla.
        grupo : array de ceros (un solo grupo, "codo").
        divisiones : dict con OD, R, angulo, N, N_filas, N_filas_L.
    """
    nodos, elementos, N_filas, N_filas_L = _malla_codo(OD, R, angulo, N)
    grupo = np.zeros(len(elementos), dtype=int)

    divisiones = dict(OD=OD, R=R, angulo=angulo, N=N, N_filas=N_filas, N_filas_L=N_filas_L)

    if verbose:
        print(
            f"Codo: OD={OD}, R={R}, angulo={angulo}°  ->  {len(nodos)} nodos, "
            f"{len(elementos)} elementos"
        )

    return dict(nodos=nodos, elementos=elementos, grupo=grupo, divisiones=divisiones)


if __name__ == "__main__":
    resultado = generar_malla_codo(OD=200, R=300, angulo=90, N=25)
    print(resultado["divisiones"])
