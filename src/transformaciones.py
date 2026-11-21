"""
transformaciones.py — helpers geometricos para posicionar piezas en 3D.

rotar_eje() es la rotación general (fórmula de Rodrigues) que usa red.py
para calcular hacia dónde apunta cada tramo (a partir de Alpha/Angulo/Giro).
traslacion()/escalado() son los mismos helpers que estaban en funciones.py,
para cuando el ensamblador tenga que mover una malla ya generada a su
posición final en la red.
"""

import numpy as np


def rotar_eje(v, eje, angulo_grados):
    """Rota el vector (o array de vectores) `v` alrededor del eje unitario
    `eje`, por `angulo_grados` grados, con la fórmula de rotación de
    Rodrigues. `eje` no necesita venir normalizado.

    La fórmula de Rodrigues descompone v en una parte paralela al eje
    (que no cambia al rotar) y una parte perpendicular (que sí rota):

        v_rot = v*cos(theta) + (eje x v)*sin(theta) + eje*(eje.v)*(1-cos(theta))
                 |__________________________|   |_____________________________|
                  mezcla v con su "giro" 90°        parte de v paralela al eje
                  alrededor del eje (cruz)           (queda igual, se suma tal cual)

    Es la misma idea que Giro() en la versión vieja (funciones.py), pero
    esa rotaba siempre alrededor del eje Y fijo; esta rota alrededor de
    cualquier eje, que es lo que necesitamos porque cada tramo de la red
    puede apuntar hacia cualquier lado.
    """
    v = np.asarray(v, dtype=float)
    eje = np.asarray(eje, dtype=float)
    eje = eje / np.linalg.norm(eje)  # por si lo pasan sin normalizar
    theta = np.deg2rad(angulo_grados)
    return (
        v * np.cos(theta)
        + np.cross(eje, v) * np.sin(theta)
        + eje * np.dot(eje, v) * (1 - np.cos(theta))
    )


def traslacion(nodos, vector):
    """Traslada un array de nodos (N, 3) sumando `vector` (broadcast)."""
    return nodos + vector


def escalado(nodos, factor):
    """Escala un array de nodos (N, 3). `factor` puede ser escalar o (x,y,z)."""
    return nodos * factor


def _matriz_rotacion(eje, angulo_grados):
    """Matriz de rotación 3x3 (Rodrigues) alrededor de `eje` (no necesita
    venir normalizado), por `angulo_grados` grados."""
    eje = np.asarray(eje, dtype=float)
    eje = eje / np.linalg.norm(eje)
    theta = np.deg2rad(angulo_grados)
    K = np.array(
        [
            [0, -eje[2], eje[1]],
            [eje[2], 0, -eje[0]],
            [-eje[1], eje[0], 0],
        ]
    )
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)


def alinear_pieza(nodos_local, eje_entrada_local, eje_normal_local,
                   eje_entrada_global, eje_normal_global):
    """Rota (sin trasladar) la malla local de una pieza (tubo/codo/tee/...)
    para que quede orientada en el espacio global.

    Cada pieza se genera en un sistema de coordenadas propio donde se sabe
    de antemano hacia qué eje local apunta su entrada (`eje_entrada_local`,
    p. ej. +X para tubo.py/reduccion.py/tee.py, +Z para codo.py) y, si la
    pieza se desvía en algún plano (un Codo, o la derivación de una Tee),
    cuál es el eje local perpendicular a ese plano (`eje_normal_local`,
    p. ej. +Y en ambos casos).

    red.py ya calculó hacia dónde debe apuntar la entrada en el espacio
    global (`eje_entrada_global`, el dir_in del nodo) y, para Codo/Tee,
    cuál es el eje de desvío global correspondiente (`eje_normal_global`,
    lo que devuelve red._eje_desvio()). Esta función arma la rotación que
    lleva el par (eje_entrada_local, eje_normal_local) exactamente a
    (eje_entrada_global, eje_normal_global), en dos pasos:

        1) una rotación mínima que lleva eje_entrada_local a
           eje_entrada_global (eje de rotación = producto cruz entre
           ambos).
        2) una rotación adicional alrededor de eje_entrada_global (ya
           alineado) que termina de llevar eje_normal_local a
           eje_normal_global (el "giro"/twist que faltaba).

    Para una pieza sin desvío (tubo, reducción) alcanza con el paso 1:
    pasá cualquier vector perpendicular a eje_entrada_local/global como
    eje_normal_local/global (el twist no importa, la sección es circular).

    Devuelve los nodos rotados (mismo shape que nodos_local), sin trasladar.
    """
    nodos_local = np.asarray(nodos_local, dtype=float)
    e_in_l = np.asarray(eje_entrada_local, dtype=float)
    e_in_l = e_in_l / np.linalg.norm(e_in_l)
    e_in_g = np.asarray(eje_entrada_global, dtype=float)
    e_in_g = e_in_g / np.linalg.norm(e_in_g)

    # --- paso 1: rotación mínima que lleva e_in_l a e_in_g ---
    v = np.cross(e_in_l, e_in_g)
    s = np.linalg.norm(v)
    c = np.dot(e_in_l, e_in_g)
    if s < 1e-9:
        if c > 0:
            R1 = np.eye(3)  # ya están alineados
        else:
            # opuestos (180°): cualquier eje perpendicular sirve
            perp = np.cross(e_in_l, [1.0, 0.0, 0.0])
            if np.linalg.norm(perp) < 1e-6:
                perp = np.cross(e_in_l, [0.0, 1.0, 0.0])
            R1 = _matriz_rotacion(perp, 180.0)
    else:
        eje1 = v / s
        angulo1 = np.degrees(np.arctan2(s, c))
        R1 = _matriz_rotacion(eje1, angulo1)

    # --- paso 2: twist alrededor de e_in_g para alinear el eje normal ---
    normal_tras_R1 = R1 @ np.asarray(eje_normal_local, dtype=float)
    e_norm_g = np.asarray(eje_normal_global, dtype=float)

    def _proyectar_perp(v, eje):
        v = v - np.dot(v, eje) * eje
        return v / np.linalg.norm(v)

    p1 = _proyectar_perp(normal_tras_R1, e_in_g)
    p2 = _proyectar_perp(e_norm_g, e_in_g)
    s2 = np.dot(np.cross(p1, p2), e_in_g)
    c2 = np.dot(p1, p2)
    angulo2 = np.degrees(np.arctan2(s2, c2))
    R2 = _matriz_rotacion(e_in_g, angulo2)

    R = R2 @ R1
    return nodos_local @ R.T
