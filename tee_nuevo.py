"""Malla preformada de una tee: tubo directo (eje x) + derivación a Alpha_2 grados.

Pasos (ver secciones numeradas más abajo):
  1. Parámetros
  2. Funciones auxiliares
  3. Tubo lateral (derivación): se recortan las caras que quedan dentro del directo
  4. Tubo directo: se recorta el hueco de la derivación y se rellena la franja libre
     hasta la curva de intersección con elementos trapezoidales
  5. Gráfico
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from funciones import FuncionesGeometria as f


# =============================================================================
# 1. PARÁMETROS
# =============================================================================
OD_1 = 200  # OD tubo directo
OD_2 = 180  # OD derivación (condición: OD_2 <= OD_1)
H_1 = 200  # de entrada a intersección
H_2 = 400  # de intersección a salida lateral
H_3 = 350  # de intersección a salida final
Alpha_2 = 60  # ángulo de la derivación [grados]

L_directo = H_1 + H_3  # longitud del tubo directo
L_salida = H_2  # longitud de la derivación

# divisiones circunferenciales de cada cilindro
N_directo = 100
N_salida = 50

N_PUNTOS_CURVA = 4000  # resolución de la curva de intersección (solo cálculo/gráfico)
MOSTRAR_LATERAL = False  # True para dibujar también el tubo lateral


# =============================================================================
# 2. FUNCIONES AUXILIARES
# =============================================================================
def compactar(nodos, elementos):
    """Elimina los nodos que ninguna cara usa y renumera la conectividad.

    Devuelve (nodos_nuevos, elementos_nuevos, nuevo_id), donde nuevo_id traduce
    numeración vieja -> nueva (-1 = nodo eliminado).
    """
    # índices (numeración vieja) de los nodos usados por alguna cara, sin repetir
    usados = np.unique(elementos)
    # tabla vieja -> nueva; -1 significa "este nodo se eliminó"
    nuevo_id = -np.ones(len(nodos), dtype=int)
    # numeración consecutiva solo para los que quedan
    nuevo_id[usados] = np.arange(len(usados))
    # coordenadas de los sobrevivientes y caras reescritas con la numeración nueva
    return nodos[usados], nuevo_id[elementos], nuevo_id


def curva_interseccion(od_1, od_2, alpha_deg, n=N_PUNTOS_CURVA):
    """Curva cerrada donde la derivación corta al tubo directo (lado z > 0).

    Se parametriza con el ángulo phi del cilindro lateral:
        y = r2 cos(phi),   z = sqrt(r1² - y²),   x = (z cos(a) - r2 sin(phi)) / sin(a)
    phi en (0, pi) da el lado cercano a la entrada y phi en (pi, 2 pi) el lejano.
    Devuelve un array (n, 3) con los puntos en orden.
    """
    r1, r2 = od_1 / 2, od_2 / 2
    a = np.radians(alpha_deg)
    phi = np.linspace(0, 2 * np.pi, n, endpoint=False)
    y = r2 * np.cos(phi)
    z = np.sqrt(r1**2 - y**2)
    x = (z * np.cos(a) - r2 * np.sin(phi)) / np.sin(a)
    return np.column_stack([x, y, z])


def caras_cruzadas_por_curva(curva, x_anillos, n_columnas, n_elementos):
    """Máscara de caras del tubo por las que pasa la curva.

    Tube numera la cara (anillo i, columna j) como i * N + j; la columna j cubre
    los ángulos [j, j+1] * 2pi/N.
    """
    i = np.clip(np.searchsorted(x_anillos, curva[:, 0]) - 1, 0, len(x_anillos) - 2)
    theta = np.mod(np.arctan2(curva[:, 2], curva[:, 1]), 2 * np.pi)
    j = np.minimum((theta / (2 * np.pi / n_columnas)).astype(int), n_columnas - 1)
    cruzadas = np.zeros(n_elementos, dtype=bool)
    cruzadas[i * n_columnas + j] = True
    return cruzadas


def aristas(elementos):
    """Aristas (nodo_a, nodo_b) de cuadriláteros en orden perimetral, ordenadas."""
    e = np.vstack(
        [
            elementos[:, [0, 1]],
            elementos[:, [1, 2]],
            elementos[:, [2, 3]],
            elementos[:, [3, 0]],
        ]
    )
    return np.sort(e, axis=1)


def lazo_de_borde(elementos_conservados, elementos_descartados):
    """Nodos (en orden) del borde del hueco: aristas compartidas entre una cara
    conservada y una descartada. Debe formar un único lazo cerrado."""
    descartadas = {tuple(a) for a in aristas(elementos_descartados)}
    borde = [
        tuple(a)
        for a in np.unique(aristas(elementos_conservados), axis=0)
        if tuple(a) in descartadas
    ]
    vecinos = {}
    for a, b in borde:
        vecinos.setdefault(a, []).append(b)
        vecinos.setdefault(b, []).append(a)
    if any(len(v) != 2 for v in vecinos.values()):
        raise ValueError("El borde del hueco no es un lazo simple (nodo con != 2 aristas)")

    inicio = borde[0][0]
    lazo = [inicio]
    previo, actual = None, inicio
    while True:
        a, b = vecinos[actual]
        siguiente = a if a != previo else b
        if siguiente == inicio:
            break
        lazo.append(siguiente)
        previo, actual = actual, siguiente
    if len(lazo) != len(vecinos):
        raise ValueError("El borde del hueco tiene más de un lazo")
    return np.array(lazo)


def trapecios_hasta_la_curva(nodos, lazo, curva):
    """Genera un nodo sobre la curva por cada nodo del borde del hueco.

    Cada nodo del borde se lleva al punto más cercano de la curva. Como varios nodos
    consecutivos pueden caer sobre el mismo punto (o en orden inverso), se corrigen
    los avances entre puntos: nunca negativos, con un mínimo, y escalados para que
    el lazo cierre exactamente una vuelta. Así los puntos avanzan siempre en el mismo
    sentido y los cuadriláteros no se cruzan ni degeneran en triángulos.
    Devuelve (lazo_orientado, nodos_curva).
    """
    n_lazo, n_curva = len(lazo), len(curva)
    P = nodos[lazo]
    idx = ((P[:, None, :] - curva[None, :, :]) ** 2).sum(axis=-1).argmin(axis=1)

    def pasos(i):  # avance entre puntos consecutivos, envuelto a [-n/2, n/2)
        return (np.diff(i, append=i[0]) + n_curva // 2) % n_curva - n_curva // 2

    if pasos(idx).sum() < 0:  # el lazo recorre la curva al revés: lo invierto
        lazo, idx = lazo[::-1], idx[::-1]

    # avances corregidos: >= 0, con un mínimo, y escalados para sumar una vuelta completa
    separacion = max(1, n_curva // (4 * n_lazo))
    avance = np.maximum(pasos(idx), 0) + separacion
    avance = avance * (n_curva / avance.sum())
    posicion = idx[0] + np.concatenate([[0], np.cumsum(avance)[:-1]])
    idx_curva = np.round(posicion).astype(int) % n_curva
    if len(np.unique(idx_curva)) != n_lazo:
        raise ValueError("Poca resolución de curva: aumentar N_PUNTOS_CURVA")
    return lazo, curva[idx_curva]


def sentido_radial(nodos, elementos):
    """Signo de (normal de la cara · dirección radial); sirve para unificar orientaciones."""
    p = nodos[elementos]
    normal = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
    centro = p.mean(axis=1)
    radial = np.column_stack([np.zeros(len(centro)), centro[:, 1], centro[:, 2]])
    return np.sign(np.einsum("ij,ij->i", normal, radial))


# =============================================================================
# 3. TUBO LATERAL (derivación)
# =============================================================================
# genero el tubo y lo giro para que quede en la posición correcta
nodos_s, elementos_s = f.Tube(OD_2, L_salida, N_salida)
nodos_s = f.Giro(nodos_s, np.radians(Alpha_2))
# conectividad como enteros y en orden perimetral
elementos_s = elementos_s.astype(int)[:, [0, 1, 3, 2]]

# True = el nodo se conserva (está fuera del cilindro directo)
mask_nodos_s = nodos_s[:, 1] ** 2 + nodos_s[:, 2] ** 2 > OD_1**2 / 4
# mask_nodos_s[elementos_s] da (M, 4): para cada cara, qué nodos se conservan;
# .all(axis=1) deja solo las caras con los 4 nodos conservados.
elementos_s = elementos_s[mask_nodos_s[elementos_s].all(axis=1)]
nodos_s, elementos_s, _ = compactar(nodos_s, elementos_s)


# =============================================================================
# 4. TUBO DIRECTO
# =============================================================================
nodos_d, elementos_d = f.Tube(OD_1, L_directo, N_directo)
nodos_d = f.traslacion(nodos_d, np.array([-H_1, 0, 0]))
elementos_d = elementos_d.astype(int)[:, [0, 1, 3, 2]]

# --- 4.1 nodos del directo que quedan dentro del cilindro lateral ---------------
# giro los nodos -Alpha_2 solo para evaluar la condición
nodos_dm = f.Giro(nodos_d, np.radians(-Alpha_2))
dentro = (nodos_dm[:, 1] ** 2 + nodos_dm[:, 2] ** 2 < OD_2**2 / 4) & (
    nodos_dm[:, 0] > 0  # lado de la derivación
)
mask_nodos_d = ~dentro  # True = el nodo se queda

# --- 4.2 caras a descartar --------------------------------------------------------
# Se descartan las caras con algún nodo adentro y también las que atraviesa la curva
# (en la punta del hueco puede pasar la curva sin que haya ningún nodo adentro).
curva = curva_interseccion(OD_1, OD_2, Alpha_2)
x_anillos = nodos_d[::N_directo, 0]  # x de cada anillo de nodos
cara_con_nodo_dentro = ~mask_nodos_d[elementos_d].all(axis=1)
cara_cruzada = caras_cruzadas_por_curva(
    curva, x_anillos, N_directo, len(elementos_d)
)
conservar = ~(cara_con_nodo_dentro | cara_cruzada)

lazo_viejo = lazo_de_borde(elementos_d[conservar], elementos_d[~conservar])
nodos_d, elementos_grilla, nuevo_id = compactar(nodos_d, elementos_d[conservar])
lazo = nuevo_id[lazo_viejo]  # lazo del hueco con la numeración nueva

# --- 4.3 relleno trapezoidal entre el borde del hueco y la curva ------------------
lazo, nodos_curva = trapecios_hasta_la_curva(nodos_d, lazo, curva)
n_lazo = len(lazo)
id_curva = len(nodos_d) + np.arange(n_lazo)  # ids de los nodos nuevos sobre la curva
siguiente = (np.arange(n_lazo) + 1) % n_lazo
elementos_relleno = np.column_stack(
    [lazo, lazo[siguiente], id_curva[siguiente], id_curva]
)

nodos_d = np.vstack([nodos_d, nodos_curva])
# misma orientación (normal hacia afuera) que las caras de la grilla
flip = sentido_radial(nodos_d, elementos_relleno) != np.median(
    sentido_radial(nodos_d, elementos_grilla)
)
elementos_relleno[flip] = elementos_relleno[flip][:, ::-1]

elementos_d = np.vstack([elementos_grilla, elementos_relleno])
n_grilla = len(elementos_grilla)
es_nodo_curva = np.arange(len(nodos_d)) >= len(nodos_d) - n_lazo


# =============================================================================
# 5. GRÁFICO
# =============================================================================
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

# coordenadas de cada vértice de cada cara: (N_elementos, 4, 3)
caras_d = nodos_d[elementos_d]
colores = ["lightsteelblue"] * n_grilla + ["khaki"] * len(elementos_relleno)
ax.add_collection3d(
    Poly3DCollection(caras_d, facecolors=colores, edgecolors="k", linewidths=0.6, alpha=0.9)
)
ax.scatter(*nodos_d[~es_nodo_curva].T, color="r", s=10)  # nodos de la grilla
ax.scatter(*nodos_d[es_nodo_curva].T, color="g", s=25)  # nodos sobre la curva
ax.plot(*np.vstack([curva, curva[:1]]).T, color="b", linewidth=1.5, label="Curva de intersección")

todos = nodos_d
if MOSTRAR_LATERAL:
    caras_s = nodos_s[elementos_s]
    ax.add_collection3d(
        Poly3DCollection(caras_s, facecolors="lightsteelblue", edgecolors="k", linewidths=0.6, alpha=0.9)
    )
    ax.scatter(*nodos_s.T, color="r", s=10)
    todos = np.vstack([nodos_d, nodos_s])

# límites y proporción real (x = eje del tubo directo)
mins, maxs = todos.min(axis=0), todos.max(axis=0)
ax.set_xlim(mins[0], maxs[0])
ax.set_ylim(mins[1], maxs[1])
ax.set_zlim(mins[2], maxs[2])
ax.set_box_aspect(maxs - mins)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()
