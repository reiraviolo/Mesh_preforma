"""
mermaid_reader.py — lector del diagrama Mermaid que describe la red de cañerías.

Lee un archivo .md/.mmd con un flowchart de Mermaid y devuelve la topología
cruda: los nodos (accesorios) con sus parámetros, y los tramos (edges) con
los suyos. No interpreta geometría ni genera ninguna malla — eso es trabajo
de red.py y del ensamblador.

Sintaxis esperada (ver Entada_de_prueba.md / prueba_red.md):

    1((Libre))
    2["Tee:OD2=100,Alpha=90,H1=200,H2=200,H3=200,Giro=90"]
    3["Codo:R=225,Angulo=90,Giro=0"]
    4["Reduccion:OD1=150,OD2=100,L=250"]

    1 -->|L=1000,OD=150| 2
    2 -->|L=2000,OD=150,S| 3
    3 -->|L=1000,OD=150,D| 4

- Nodo libre:        id((Tipo))                  (típicamente Tipo = "Libre")
- Nodo con params:   id["Tipo:clave=valor,..."]
- Tramo:             origen -->|clave=valor,...,flag| destino

Un parámetro de tramo/nodo sin "=" (por ejemplo "S" o "D") se guarda como
flag booleano: {"S": True}.
"""

import re

PATRON_NODO_LIBRE = re.compile(r"^\s*(\w+)\(\((\w+)\)\)\s*$")
PATRON_NODO_PARAMS = re.compile(r'^\s*(\w+)\["([^:]+):([^"]+)"\]\s*$')
PATRON_TRAMO = re.compile(r"^\s*(\w+)\s*-->\s*(?:\|([^|]*)\|)?\s*(\w+)\s*;?\s*$")

LINEAS_IGNORAR = ("```", "flowchart", "graph")


def _parsear_params(texto):
    """Convierte 'OD2=100,Alpha=90,S' en {'OD2': 100.0, 'Alpha': 90.0, 'S': True}.

    Acepta tanto 'clave=valor' como 'clave:valor' (algunos archivos usan
    ':' en vez de '='; se tratan igual para no romper por un detalle de
    puntuación)."""
    params = {}
    if not texto:
        return params
    for token in texto.split(","):
        token = token.strip()
        if not token:
            continue
        if "=" in token or ":" in token:
            clave, valor = re.split(r"[=:]", token, maxsplit=1)
            clave = clave.strip()
            valor = valor.strip()
            try:
                valor_num = float(valor)
                if valor_num == int(valor_num):
                    valor_num = int(valor_num)
                params[clave] = valor_num
            except ValueError:
                params[clave] = valor  # queda como texto si no es numérico
        else:
            params[token] = True
    return params


def parsear_mermaid(texto):
    """
    Parsea el texto de un Mermaid (ya en memoria, no un archivo) y devuelve
    (nodos, tramos). Es lo que usa leer_mermaid() después de leer el
    archivo; separado así para poder parsear directo el contenido de un
    editor de texto (p.ej. la app gráfica) sin pasar por disco.

    nodos : dict
        {id: {"tipo": str, "params": dict}}
    tramos : list
        [{"origen": id, "destino": id, "params": dict}, ...]
    """
    nodos = {}
    tramos = []

    for linea in texto.splitlines():
        contenido = linea.strip()
        if not contenido:
            continue
        if any(contenido.startswith(p) for p in LINEAS_IGNORAR):
            continue

        m = PATRON_NODO_LIBRE.match(contenido)
        if m:
            nodo_id, tipo = m.groups()
            nodos[nodo_id] = {"tipo": tipo, "params": {}}
            continue

        m = PATRON_NODO_PARAMS.match(contenido)
        if m:
            nodo_id, tipo, params_txt = m.groups()
            nodos[nodo_id] = {"tipo": tipo.strip(), "params": _parsear_params(params_txt)}
            continue

        m = PATRON_TRAMO.match(contenido)
        if m:
            origen, params_txt, destino = m.groups()
            tramos.append(
                {
                    "origen": origen,
                    "destino": destino,
                    "params": _parsear_params(params_txt),
                }
            )
            continue

        # Línea que no matcheó ningún patrón conocido: la ignoramos pero
        # avisamos, para no tragarnos silenciosamente un error de sintaxis.
        print(f"[mermaid_reader] línea no reconocida, ignorada: {contenido!r}")

    # Cualquier nodo mencionado en un tramo pero sin línea de definición
    # propia se agrega como "Libre" por defecto.
    for tramo in tramos:
        for extremo in (tramo["origen"], tramo["destino"]):
            if extremo not in nodos:
                nodos[extremo] = {"tipo": "Libre", "params": {}}
                print(f"[mermaid_reader] nodo {extremo!r} sin definición propia, asumido Libre")

    return nodos, tramos


def leer_mermaid(path):
    """Lee un archivo Mermaid de disco y lo parsea (ver parsear_mermaid)."""
    with open(path, "r", encoding="utf-8") as f:
        return parsear_mermaid(f.read())


if __name__ == "__main__":
    import sys
    import pprint

    if len(sys.argv) < 2:
        print("Uso: python mermaid_reader.py archivo.md")
        sys.exit(1)

    nodos, tramos = leer_mermaid(sys.argv[1])
    print("NODOS:")
    pprint.pprint(nodos)
    print("\nTRAMOS:")
    pprint.pprint(tramos)
