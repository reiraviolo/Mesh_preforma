"""
main.py — punto de entrada de la app. Lee un archivo Mermaid con la
topología de la red de cañerías, genera cada pieza (tubo/codo/reducción/
tee) ya ubicada en el espacio, grafica el resultado, funde las mallas en
una sola y exporta el archivo .malla final.

Uso:
    uv run python main.py archivo.md [N] [salida.malla]

    archivo.md   : ruta al Mermaid (ver prueba_red.md / Caso_1.md como
                   ejemplo de sintaxis).
    N            : divisiones de perímetro del tronco principal (por
                   defecto 24). En la derivación de cada Tee, tee.py
                   elige su propio N y ese es el que usan el tramo y los
                   accesorios siguientes de esa rama (ver ensamblador.py).
    salida.malla : ruta del archivo .malla a generar (por defecto, mismo
                   nombre que archivo.md pero con extensión .malla).
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, "src")
sys.path.insert(0, "src/parts")

from mermaid_reader import leer_mermaid
from red import calcular_posiciones
from ensamblador import ensamblar, fusionar_piezas
from exportador import escribir_malla
import graficador
import tee


def main():
    if len(sys.argv) < 2:
        print("Uso: uv run python main.py archivo.md [N] [salida.malla]")
        sys.exit(1)

    archivo = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    if len(sys.argv) > 3:
        salida_malla = sys.argv[3]
    else:
        salida_malla = os.path.splitext(archivo)[0] + ".malla"

    nodos, tramos = leer_mermaid(archivo)
    piezas = ensamblar(nodos, tramos, N=N)
    posiciones = calcular_posiciones(nodos, tramos)

    # --- fusion de piezas en una sola malla (nodos coincidentes unificados) ---
    # zonas_soldadura: lista de arrays de ids de nodo, uno por anillo donde
    # dos piezas se fusionaron (= donde, en la caneria real, va la soldadura)
    nodos_f, elementos_f, zonas_soldadura = fusionar_piezas(piezas)
    control = tee.control_malla(nodos_f, elementos_f)
    print("Control de malla fusionada:", control)
    print(f"Zonas de soldadura: {len(zonas_soldadura)} "
          f"({sum(len(z) for z in zonas_soldadura)} nodos de union)")

    mermaid_texto = open(archivo, encoding="utf-8").read()
    escribir_malla(salida_malla, nodos_f, elementos_f, mermaid_texto=mermaid_texto,
                   zonas_soldadura=zonas_soldadura)
    print(f"Malla exportada: {salida_malla} "
          f"({len(nodos_f)} nodos, {len(elementos_f)} elementos)")

    # Se grafica la malla YA FUSIONADA (un solo conjunto de nodos/elementos),
    # no cada pieza por separado: asi en las uniones (tee/codo/reduccion) no
    # queda la doble linea que se veia al dibujar cada pieza con su propio
    # Poly3DCollection (dos superficies casi-coincidentes superpuestas).
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    graficador.graficar_malla(nodos_f, elementos_f, ax=ax, mostrar=False,
                               colores={0: "#ffcf80"})

    # marca en rojo los nodos de soldadura (uniones entre piezas)
    graficador.graficar_soldaduras(nodos_f, zonas_soldadura, ax=ax, mostrar=False)

    # superpone las lineas de eje/centro de cada tramo sobre la malla
    graficador.graficar_ejes(posiciones, tramos, nodos, ax=ax, mostrar=False)

    # zoom al conjunto completo (graficar_malla ya ajusta los limites con
    # nodos_f, que incluye todas las piezas, asi que alcanza con un 5% de
    # margen extra para que no quede justo al borde)
    lo, hi = nodos_f.min(axis=0), nodos_f.max(axis=0)
    centro, radio = (lo + hi) / 2, max((hi - lo).max() / 2, 1e-9) * 1.05
    ax.set_xlim(centro[0] - radio, centro[0] + radio)
    ax.set_ylim(centro[1] - radio, centro[1] + radio)
    ax.set_zlim(centro[2] - radio, centro[2] + radio)
    ax.set_box_aspect((1, 1, 1))

    ax.set_title(f"{archivo} — {len(piezas)} piezas, N={N} — malla unificada")
    plt.show()


if __name__ == "__main__":
    main()
