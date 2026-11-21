"""
exportador.py — guarda una malla (nodos + elementos) en un archivo de
texto plano ".malla", con el Mermaid de origen como encabezado y, además,
las zonas de soldadura: los nodos de unión entre piezas (donde dos piezas
-tubo/codo/reducción/tee- se fusionan y que, en la cañería real, son la
soldadura entre ambas).

Formato del archivo:

    ### MERMAID ###
    # <cada línea del mermaid de origen, comentada con "# ">
    ### NODOS ###
    N <cantidad de nodos>
    <id> <x> <y> <z>
    ...
    ### ELEMENTOS ###
    E <cantidad de elementos>
    <n1> <n2> <n3> <n4>
    ...
    ### SOLDADURAS ###
    Z <cantidad de zonas de soldadura>
    <zona> <cantidad de nodos de la zona> <id1> <id2> ... <idN>
    ...

Un triángulo se guarda como cuadrilátero con el último nodo repetido
(convención ya usada en todo el proyecto: [a, b, c, c]).

Cada "zona" de soldadura es un anillo de nodos: la circunferencia donde
una pieza queda soldada a la siguiente. La sección ### SOLDADURAS ### es
opcional al leer (archivos .malla viejos, sin esa sección, se siguen
leyendo bien y devuelven zonas_soldadura=[]).
"""

import numpy as np


def escribir_malla(path, nodos, elementos, mermaid_texto="", zonas_soldadura=None):
    """
    Escribe la malla fusionada a `path` (típicamente con extensión .malla).

    nodos            : (N, 3) float
    elementos        : (E, 4) int
    mermaid_texto    : el contenido del archivo Mermaid que generó esta red,
                       se guarda tal cual como comentario de encabezado.
    zonas_soldadura  : lista de arrays de ids de nodo (uno por zona/anillo
                       de soldadura). Si es None, no se escribe la sección.
    """
    nodos = np.asarray(nodos, dtype=float)
    elementos = np.asarray(elementos, dtype=int)

    with open(path, "w", encoding="utf-8") as f:
        f.write("### MERMAID ###\n")
        for linea in mermaid_texto.splitlines():
            f.write(f"# {linea}\n")

        f.write("### NODOS ###\n")
        f.write(f"N {len(nodos)}\n")
        for i, (x, y, z) in enumerate(nodos):
            f.write(f"{i} {x:.6f} {y:.6f} {z:.6f}\n")

        f.write("### ELEMENTOS ###\n")
        f.write(f"E {len(elementos)}\n")
        for e in elementos:
            f.write(" ".join(str(int(n)) for n in e) + "\n")

        if zonas_soldadura is not None:
            f.write("### SOLDADURAS ###\n")
            f.write(f"Z {len(zonas_soldadura)}\n")
            for z, nodos_zona in enumerate(zonas_soldadura):
                ids = " ".join(str(int(n)) for n in nodos_zona)
                f.write(f"{z} {len(nodos_zona)} {ids}\n")


def leer_malla(path):
    """Lee un archivo .malla escrito con escribir_malla(). Devuelve
    (nodos, elementos, mermaid_texto, zonas_soldadura).

    zonas_soldadura es una lista de arrays de ids de nodo (uno por zona);
    si el archivo no tiene la sección ### SOLDADURAS ### (formato viejo),
    se devuelve [].
    """
    with open(path, "r", encoding="utf-8") as f:
        lineas = [l.rstrip("\n") for l in f]

    i = 0
    assert lineas[i].strip() == "### MERMAID ###"
    i += 1
    mermaid_lineas = []
    while lineas[i].strip() != "### NODOS ###":
        mermaid_lineas.append(lineas[i][2:] if lineas[i].startswith("# ") else lineas[i])
        i += 1
    mermaid_texto = "\n".join(mermaid_lineas)

    i += 1  # saltar "### NODOS ###"
    n_nodos = int(lineas[i].split()[1])
    i += 1
    nodos = np.zeros((n_nodos, 3))
    for k in range(n_nodos):
        partes = lineas[i].split()
        nodos[k] = [float(partes[1]), float(partes[2]), float(partes[3])]
        i += 1

    assert lineas[i].strip() == "### ELEMENTOS ###"
    i += 1
    n_elem = int(lineas[i].split()[1])
    i += 1
    elementos = np.zeros((n_elem, 4), dtype=int)
    for k in range(n_elem):
        elementos[k] = [int(x) for x in lineas[i].split()]
        i += 1

    zonas_soldadura = []
    if i < len(lineas) and lineas[i].strip() == "### SOLDADURAS ###":
        i += 1
        n_zonas = int(lineas[i].split()[1])
        i += 1
        for _ in range(n_zonas):
            partes = lineas[i].split()
            # partes[0]=zona, partes[1]=cantidad, el resto son los ids
            zonas_soldadura.append(np.array([int(x) for x in partes[2:]], dtype=int))
            i += 1

    return nodos, elementos, mermaid_texto, zonas_soldadura
