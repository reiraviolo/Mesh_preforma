"""Lectura del archivo de entrada (.pmsh).

Ejemplo de archivo:

    ###PARAMETROS######
    N=25

    ###TRAMOS######
    TRAMO;DIR(X,Y,Z);LONG;OD;Nodo_I;Nodo_F
    1;1,0,0;1000;150;1;2

    ###NODOS######
    1;;
    2;Tee;150;150;90;200;200;200
"""

import sys
import numpy as np
from modelo import Nodo, Red, Tramo

# Nombre de cada parámetro de los accesorios, en el orden en que aparecen en el txt
CAMPOS = {
    "Tee": ["OD1", "OD2", "alpha", "H1", "H2", "H3"],
    "Codo": ["OD", "angulo", "R", "H"],
    "Reduccion": ["OD1", "OD2", "L"],
}


def leer_red(ruta):
    """Lee el archivo y devuelve un objeto Red con los tramos y los nodos."""
    params = {"N": 25}  # valor por defecto si el archivo no lo define
    tramos = {}
    nodos = {}

    seccion = None  # sección en la que estoy parado: "PARAMETROS", "TRAMOS" o "NODOS"

    with open(ruta, encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()  # saco espacios y el salto de línea

            # Línea vacía: la salteo
            if linea == "":
                continue

            # Línea de sección, por ejemplo "###TRAMOS######"
            if linea.startswith("###"):
                if "PARAMETROS" in linea:
                    seccion = "PARAMETROS"
                elif "TRAMOS" in linea:
                    seccion = "TRAMOS"
                elif "NODOS" in linea:
                    seccion = "NODOS"
                continue

            # Línea de datos: depende de la sección
            if seccion == "PARAMETROS":
                clave, valor = linea.split("=")
                params[clave.strip()] = float(valor)

            elif seccion == "TRAMOS":
                if linea.startswith("TRAMO"):  # encabezado de columnas
                    continue
                tramo = leer_tramo(linea)
                tramos[tramo.id] = tramo

            elif seccion == "NODOS":
                nodo = leer_nodo(linea)
                nodos[nodo.id] = nodo

    params["N"] = int(params["N"])  # N es una cantidad de divisiones: entero
    return Red(tramos=tramos, nodos=nodos, params=params)


def leer_tramo(linea):
    """'1;1,0,0;1000;150;1;2' -> Tramo"""
    campos = linea.split(";")

    direccion = np.array(campos[1].split(","), dtype=float)
    direccion = direccion / np.linalg.norm(direccion)  # la normalizo (largo 1)

    return Tramo(
        id=int(campos[0]),
        dir=direccion,
        long=float(campos[2]),
        od=float(campos[3]),
        nodo_i=int(campos[4]),
        nodo_f=int(campos[5]),
    )


def leer_nodo(linea):
    """'2;Tee;150;150;90;200;200;200' -> Nodo"""
    campos = linea.split(";")
    id_nodo = int(campos[0])
    tipo = campos[1]

    # Nodo sin accesorio (extremo libre), por ejemplo '1;;'
    if tipo == "":
        return Nodo(id=id_nodo, tipo=None)

    # Asocio cada valor con su nombre según CAMPOS
    nombres = CAMPOS[tipo]
    valores = campos[2:]
    params = {}
    for nombre, valor in zip(nombres, valores):
        params[nombre] = float(valor)

    return Nodo(id=id_nodo, tipo=tipo, params=params)


if __name__ == "__main__":
    # Uso: python entrada.py archivo.pmsh
    red = leer_red(sys.argv[1])

    print("Parámetros:", red.params)
    for tramo in red.tramos.values():
        print(tramo)
    for nodo in red.nodos.values():
        print(nodo)
