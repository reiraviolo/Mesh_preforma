"""Código que llama a cada uno de los accesorios y genera 
el gráfico 3D de cada uno según sean llamados."""

import matplotlib.pyplot as plt

from entrada import leer_red
from funciones import FuncionesGeometria
from modelo import Nodo
from graficar import graficar


red = leer_red("ejemplo_entrada.pmsh") #leearchivo de entrada y genera la red
N = red.params["N"]

for nodo in red.nodos.values():
    p = nodo.params
    if nodo.tipo == "Codo":
        nodos, elementos = FuncionesGeometria.Codo(p["OD"], p["R"], p["angulo"], N)
        graficar(nodos, elementos, "Codo")
    elif nodo.tipo == "Reduccion":
        nodos, elementos = FuncionesGeometria.Reduccion(p["OD1"], p["OD2"], p["L"], N)
        graficar(nodos, elementos, "Reduccion")
    elif nodo.tipo == "Tee":
        ...  # todavía no hay función


plt.show()