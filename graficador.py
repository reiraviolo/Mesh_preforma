"""
graficador.py — Dibuja una malla de cuadriláteros a partir de sus matrices de nodos y elementos.

No conoce nada de la geometría de la tee: solo recibe arrays de nodos (N,3) y elementos (E,4)
(más, opcionalmente, un array de grupo para colorear). Pensado para usarse con lo que devuelve
`tee.generar_malla_tee`, pero sirve para cualquier malla de cuadriláteros con esta forma.

Un elemento triangular se admite como cuadrilátero con el último nodo repetido: [a, b, c, c].

Uso típico
----------
    from tee import generar_malla_tee
    from graficador import graficar_malla, graficar_desarrollado

    r = generar_malla_tee(OD_1=200, OD_2=100, H_1=350, H_2=350, H_3=400, Alpha_2=30)
    graficar_malla(r["nodos"], r["elementos"], r["grupo"])
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.collections import PolyCollection


COLORES_POR_GRUPO = {0: '#ffe0b3', 1: '#cfe8c8', 2: '#cfe0ff'}   # 0=derivación, 1=corona, 2=directo (convención de tee.py)
NOMBRES_GRUPO = {0: 'derivación', 1: 'corona', 2: 'directo'}


def graficar_malla(nodos, elementos, grupo=None, colores=None, edgecolor='k', linewidth=0.6, alpha=0.95,
                    elev=22, azim=-60, figsize=(10, 8), titulo=None, ax=None, mostrar=True):
    """Dibuja la malla en 3D (superficie de cuadriláteros).

    nodos     : (N,3) float
    elementos : (E,4) int — conectividad; un triángulo se pasa como [a,b,c,c]
    grupo     : (E,) int opcional, para colorear cada elemento según su grupo
    colores   : dict {valor_de_grupo: color}; por defecto usa la convención 0/1/2 de tee.py.
                Si `grupo` es None, se usa un solo color (el de la clave 0, o 'tab:blue').
    ax        : Axes3D existente, para dibujar varias mallas juntas (p. ej. media pieza + espejo)
    mostrar   : si True llama a plt.show() al final (ponelo en False si vas a seguir agregando cosas al ax)

    Devuelve (fig, ax).
    """
    nodos = np.asarray(nodos, float)
    elementos = np.asarray(elementos, int)
    colores = colores or COLORES_POR_GRUPO

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig = ax.figure

    if grupo is None:
        fc = colores.get(0, 'tab:blue')
    else:
        grupo = np.asarray(grupo, int)
        fc = [colores.get(g, 'lightgray') for g in grupo]

    ax.add_collection3d(Poly3DCollection(nodos[elementos], facecolors=fc, edgecolors=edgecolor,
                                          linewidths=linewidth, alpha=alpha))

    lo, hi = nodos.min(axis=0), nodos.max(axis=0)
    centro, radio = (lo + hi) / 2, max((hi - lo).max() / 2, 1e-9)
    ax.set_xlim(centro[0] - radio, centro[0] + radio)
    ax.set_ylim(centro[1] - radio, centro[1] + radio)
    ax.set_zlim(centro[2] - radio, centro[2] + radio)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.view_init(elev=elev, azim=azim)
    if titulo:
        ax.set_title(titulo)

    if grupo is not None:
        for g in sorted(set(grupo.tolist())):
            ax.plot([], [], 's', color=colores.get(g, 'lightgray'), label=NOMBRES_GRUPO.get(g, f'grupo {g}'))
        ax.legend(loc='upper left')

    if mostrar:
        plt.show()
    return fig, ax


def graficar_desarrollado(nodos, elementos, grupo=None, radio=None, colores=None, edgecolor='k', linewidth=0.5,
                           figsize=(14, 6), titulo=None, ax=None, mostrar=True):
    """Dibuja la malla "desarrollada" en 2D: (x, radio·atan2(z, y)) — útil para ver de un vistazo la
    regularidad de la malla sobre un caño cilíndrico de eje X y radio `radio`.

    Si `radio` es None, se estima como el radio medio de los nodos respecto del eje X (sirve como
    vista aproximada incluso si la malla no es exactamente un cilindro de eje X, p. ej. incluye la
    derivación en otro eje: esos elementos van a salir deformados en este gráfico, pero no superpuestos
    de forma engañosa si se grafican por separado).

    Devuelve (fig, ax).
    """
    nodos = np.asarray(nodos, float)
    elementos = np.asarray(elementos, int)
    colores = colores or COLORES_POR_GRUPO

    if radio is None:
        radio = float(np.hypot(nodos[:, 1], nodos[:, 2]).mean())

    uv = np.column_stack([nodos[:, 0], radio * np.arctan2(nodos[:, 2], nodos[:, 1])])

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    if grupo is None:
        ax.add_collection(PolyCollection(uv[elementos], facecolors=colores.get(0, 'tab:blue'),
                                          edgecolors=edgecolor, linewidths=linewidth))
    else:
        grupo = np.asarray(grupo, int)
        for g in sorted(set(grupo.tolist())):
            m = grupo == g
            ax.add_collection(PolyCollection(uv[elementos[m]], facecolors=colores.get(g, 'lightgray'),
                                              edgecolors=edgecolor, linewidths=linewidth,
                                              label=NOMBRES_GRUPO.get(g, f'grupo {g}')))
        ax.legend(loc='upper right')

    ax.autoscale_view()
    ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel(f'{radio:.1f} · atan2(z, y)')
    if titulo:
        ax.set_title(titulo)
    if mostrar:
        plt.show()
    return fig, ax


if __name__ == "__main__":
    # Ejemplo de uso end-to-end: genera la malla con tee.py y la dibuja.
    from tee import generar_malla_tee

    r = generar_malla_tee(OD_1=200, OD_2=100, H_1=350, H_2=350, H_3=400, Alpha_2=30, verbose=True)
    graficar_malla(r["nodos"], r["elementos"], r["grupo"],
                   titulo=f"Tee OD_1=200, OD_2=100, Alpha_2=30° — {len(r['elementos'])} elementos")
