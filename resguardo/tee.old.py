"""Malla preformada de una tee: tubo directo (eje x) + derivación a Alpha_2 grados.

Pasos (ver secciones numeradas más abajo):
  1. Parámetros
  2. Funciones auxiliares
  3. Curva de intersección entre los dos cilindros
  4. Tubo directo recortado contra la derivación
  5. Derivación: se prueban distintas divisiones circunferenciales y se queda la que
     da los elementos de la unión más cuadrados. Los dos tubos comparten los mismos
     nodos sobre la curva (coincidentes) y se unen a ella con una franja de trapecios.
  6. Gráfico
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

# divisiones circunferenciales
N_directo = 50
N_salida = 50  # valor inicial de la derivación (se ajusta si AJUSTAR_DIVISIONES)

# La derivación prueba N_salida - BUSQUEDA ... N_salida + BUSQUEDA (solo valores pares)
# y se queda con la que da los elementos de la unión más cuadrados.
AJUSTAR_DIVISIONES = True
BUSQUEDA = 12

# arista mínima admitida en el relleno = FACTOR * tamaño de elemento de la grilla
FACTOR_ARISTA_MIN = 0.5

N_PUNTOS_CURVA = 4000  # resolución de la curva de intersección (solo cálculo/gráfico)
MOSTRAR_LATERAL = True  # False para dibujar solo el tubo directo


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

    La curva debe estar expresada en el sistema propio del tubo (eje del tubo = x).
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
    """Aristas (nodo_a, nodo_b) de cuadriláteros en orden perimetral, ordenadas.

    Apila las 4 aristas de todas las caras por bloques: la arista k pertenece a la
    cara k % M.
    """
    e = np.vstack(
        [
            elementos[:, [0, 1]],
            elementos[:, [1, 2]],
            elementos[:, [2, 3]],
            elementos[:, [3, 0]],
        ]
    )
    return np.sort(e, axis=1)


def tamano_de_elemento(nodos, elementos):
    """Tamaño típico de elemento: mediana de la longitud de sus aristas."""
    a = np.unique(aristas(elementos), axis=0)
    return np.median(np.linalg.norm(nodos[a[:, 0]] - nodos[a[:, 1]], axis=1))


def distancia_a_curva(puntos, curva, bloque=500):
    """Distancia de cada punto al punto más cercano de la curva (por bloques, para
    no armar una matriz gigante cuando la malla es fina)."""
    d = np.empty(len(puntos))
    for k in range(0, len(puntos), bloque):
        p = puntos[k : k + bloque]
        d[k : k + bloque] = np.sqrt(
            ((p[:, None, :] - curva[None, :, :]) ** 2).sum(axis=-1).min(axis=1)
        )
    return d


def aristas_de_borde(elementos_conservados, elementos_descartados):
    """Aristas compartidas entre una cara conservada y una descartada (borde del hueco)."""
    descartadas = {tuple(a) for a in aristas(elementos_descartados)}
    return [
        tuple(a)
        for a in np.unique(aristas(elementos_conservados), axis=0)
        if tuple(a) in descartadas
    ]


def limpiar_hueco(elementos, conservar):
    """Deja el borde del hueco como un único lazo simple.

    Al quitar más nodos pueden quedar caras sueltas (islas unidas al resto solo por
    un nodo) o puntos donde el borde se toca a sí mismo. Se repite hasta que no
    quede ninguno:
      1. se descarta todo lo que no sea la mayor porción conectada por aristas;
      2. se descartan las caras conservadas que tocan un nodo del borde con más de
         dos aristas de borde (el borde se pellizca ahí).
    """
    conservar = conservar.copy()
    while True:
        cambio = False

        # 1) quedarse con la porción conexa más grande (caras que comparten arista)
        idx = np.flatnonzero(conservar)
        m = len(idx)
        padre = np.arange(m)

        def raiz(a):
            while padre[a] != a:
                padre[a] = padre[padre[a]]
                a = padre[a]
            return a

        primera = {}
        for k, arista in enumerate(map(tuple, aristas(elementos[idx]))):
            cara = k % m
            if arista in primera:
                padre[raiz(cara)] = raiz(primera[arista])
            else:
                primera[arista] = cara
        raices = np.array([raiz(c) for c in range(m)])
        mayor = np.bincount(raices).argmax()
        if (raices != mayor).any():
            conservar[idx[raices != mayor]] = False
            cambio = True

        # 2) nodos del borde con más de 2 aristas de borde
        borde = aristas_de_borde(elementos[conservar], elementos[~conservar])
        grado = {}
        for a, b in borde:
            grado[a] = grado.get(a, 0) + 1
            grado[b] = grado.get(b, 0) + 1
        pellizcos = [n for n, g in grado.items() if g > 2]
        if pellizcos:
            conservar[np.isin(elementos, pellizcos).any(axis=1)] = False
            cambio = True

        if not cambio:
            return conservar


def lazo_de_borde(elementos_conservados, elementos_descartados):
    """Nodos (en orden) del borde del hueco: aristas compartidas entre una cara
    conservada y una descartada. Debe formar un único lazo cerrado."""
    borde = aristas_de_borde(elementos_conservados, elementos_descartados)
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


def preparar_recorte(nodos, elementos, dentro, cara_cruzada, curva):
    """Decide qué caras de un tubo se conservan contra el otro cilindro.

    nodos, elementos : malla completa del tubo (elementos en orden perimetral)
    dentro           : máscara de nodos que están dentro del otro cilindro
    cara_cruzada     : máscara de caras por las que pasa la curva

    Se quitan también los nodos a menos de arista_min de la curva: así ninguna unión
    nodo -> curva puede ser más corta que eso. Se descartan las caras con algún nodo
    quitado o cruzadas por la curva, y el hueco queda como un único lazo simple.
    Devuelve un dict con la malla, la máscara 'conservar', 'tam' y 'arista_min'.
    """
    tam = tamano_de_elemento(nodos, elementos)
    arista_min = FACTOR_ARISTA_MIN * tam
    nodo_se_queda = ~(dentro | (distancia_a_curva(nodos, curva) < arista_min))
    cara_con_nodo_quitado = ~nodo_se_queda[elementos].all(axis=1)
    conservar = limpiar_hueco(elementos, ~(cara_con_nodo_quitado | cara_cruzada))
    return dict(
        nodos=nodos, elementos=elementos, conservar=conservar, tam=tam, arista_min=arista_min
    )


def cantidad_de_nodos_del_lazo(r):
    el, c = r["elementos"], r["conservar"]
    return len(lazo_de_borde(el[c], el[~c]))


def ampliar_lazo(r, curva, n_objetivo):
    """Agrega nodos al lazo del hueco (de a 2) descartando caras del borde.

    Una cara conservada con una sola arista sobre el borde, al descartarla, cambia esa
    arista por tres: el lazo gana 2 nodos. Se prueban primero las caras más cercanas a
    la curva (el hueco crece lo menos posible) y se acepta la primera que deja el
    borde como un único lazo.
    """
    el = r["elementos"]
    conservar = r["conservar"].copy()
    centros = r["nodos"][el].mean(axis=1)
    n = len(lazo_de_borde(el[conservar], el[~conservar]))
    if (n_objetivo - n) % 2:
        raise ValueError("La diferencia de nodos entre los lazos es impar")
    while n < n_objetivo:
        borde = set(aristas_de_borde(el[conservar], el[~conservar]))
        en_borde = np.array([tuple(a) in borde for a in aristas(el)]).reshape(4, -1).sum(axis=0)
        candidatas = np.flatnonzero(conservar & (en_borde == 1))
        candidatas = candidatas[np.argsort(distancia_a_curva(centros[candidatas], curva))]
        for c in candidatas:
            prueba = conservar.copy()
            prueba[c] = False
            try:
                m = len(lazo_de_borde(el[prueba], el[~prueba]))
            except ValueError:
                continue
            if m == n + 2:
                conservar, n = prueba, m
                break
        else:
            raise ValueError("No se pudo ampliar el lazo del hueco")
    return conservar


def lazo_final(r):
    """Malla conservada (renumerada) y su lazo de borde en la numeración nueva."""
    el, c = r["elementos"], r["conservar"]
    lazo_viejo = lazo_de_borde(el[c], el[~c])
    nodos, elementos, nuevo_id = compactar(r["nodos"], el[c])
    return nodos, elementos, nuevo_id[lazo_viejo]


def forma_cuadrilateros(p0, p1, p2, p3):
    """4 * área / suma de aristas² de cuadriláteros (p0, p1, p2, p3) en orden perimetral.

    Vale 1 para un cuadrado y baja si el elemento es alargado, achatado o torcido.
    Acepta arrays (..., 3).
    """
    area = 0.5 * np.linalg.norm(np.cross(p2 - p0, p3 - p1), axis=-1)
    suma = sum(((b - a) ** 2).sum(axis=-1) for a, b in ((p0, p1), (p1, p2), (p2, p3), (p3, p0)))
    return 4 * area / suma


def optimizar_nodos_curva(idx, curva, lazos, arista_min, barridos=8, n_cand=41):
    """Mueve cada nodo común sobre la curva (entre sus vecinos) para mejorar la forma
    de los cuadriláteros que lo rodean.

    idx         : índice en la curva de cada nodo común (en orden)
    lazos       : [P_a, P_b] posiciones de los nodos de borde de cada tubo, (n, 3) cada una
    Cada nodo toca 4 cuadriláteros (2 por tubo). Se prueban n_cand posiciones entre
    sus vecinos, respetando la distancia mínima a ellos, y se elige la que maximiza
    la peor forma (con un pequeño peso para el promedio). Se repite por barridos.
    """
    n, nc = len(idx), len(curva)
    idx = idx.copy()
    for _ in range(barridos):
        mejoro = False
        for i in range(n):
            ant, sig = idx[i - 1], idx[(i + 1) % n]
            span = (sig - ant) % nc
            if span < 3:
                continue
            cand = np.unique((ant + np.linspace(1, span - 1, n_cand).astype(int)) % nc)
            cand = np.union1d(cand, [idx[i]])
            Q = curva[cand]  # (K, 3)
            Cs, Ca = curva[idx[(i + 1) % n]], curva[idx[i - 1]]
            valido = (np.linalg.norm(Q - Cs, axis=1) >= arista_min) & (
                np.linalg.norm(Q - Ca, axis=1) >= arista_min
            )
            q = []
            for P in lazos:  # cuadriláteros i-1 e i de cada tubo
                q.append(forma_cuadrilateros(P[i - 1], P[i], Q, Ca))
                q.append(forma_cuadrilateros(P[i], P[(i + 1) % n], Cs, Q))
            q = np.array(q)  # (4, K)
            puntaje = np.where(valido, q.min(axis=0) + 0.1 * q.mean(axis=0), -np.inf)
            mejor = cand[np.argmax(puntaje)]
            actual = puntaje[np.flatnonzero(cand == idx[i])[0]]
            if puntaje.max() > actual + 1e-9 and mejor != idx[i]:
                idx[i] = mejor
                mejoro = True
        if not mejoro:
            break
    return idx


def nodos_comunes_en_la_curva(nodos_a, lazo_a, nodos_b, lazo_b, curva, arista_min):
    """Un único juego de nodos sobre la curva para los dos tubos.

    Los dos lazos deben tener la misma cantidad de nodos. Cada nodo de cada lazo se
    lleva al punto más cercano de la curva (posición medida como longitud de arco);
    los dos lazos se orientan en el mismo sentido y se alinean cíclicamente para que
    el nodo i de uno quede lo más cerca posible del nodo i del otro. La posición final
    de cada nodo común es el promedio, corregida para que:
      - los puntos avancen siempre en el mismo sentido (sin cruces ni triángulos),
      - cierren exactamente una vuelta,
      - la distancia entre nodos consecutivos sea >= arista_min (se mezcla con el
        reparto uniforme lo mínimo necesario).
    Después se optimiza la posición de cada nodo (optimizar_nodos_curva).
    Devuelve (lazo_a, lazo_b, nodos_curva, cumple) con los lazos ya orientados y
    alineados; cumple indica si todas las aristas sobre la curva son >= arista_min.
    """
    n = len(lazo_a)
    tramo = np.linalg.norm(np.roll(curva, -1, axis=0) - curva, axis=1)
    arco = np.concatenate([[0], np.cumsum(tramo)])
    perimetro = arco[-1]
    arco = arco[:-1]

    def envolver(x):  # a [-perimetro/2, perimetro/2)
        return (x + perimetro / 2) % perimetro - perimetro / 2

    def posiciones(nodos, lazo):
        P = nodos[lazo]
        idx = ((P[:, None, :] - curva[None, :, :]) ** 2).sum(axis=-1).argmin(axis=1)
        a = arco[idx]
        if envolver(np.diff(a, append=a[0])).sum() < 0:  # recorre la curva al revés
            lazo, a = lazo[::-1], a[::-1]
        return lazo, a

    lazo_a, a_a = posiciones(nodos_a, lazo_a)
    lazo_b, a_b = posiciones(nodos_b, lazo_b)

    # alineo cíclicamente b con a (el nodo i de b queda cerca del nodo i de a)
    costos = [np.sum(envolver(np.roll(a_b, -k) - a_a) ** 2) for k in range(n)]
    k = int(np.argmin(costos))
    lazo_b, a_b = np.roll(lazo_b, -k), np.roll(a_b, -k)

    objetivo = a_a + 0.5 * envolver(a_b - a_a)  # promedio de lo que pide cada tubo
    avance = np.maximum(envolver(np.diff(objetivo, append=objetivo[0])), 0)
    avance = avance + 0.05 * perimetro / n
    avance *= perimetro / avance.sum()
    uniforme = perimetro / n

    for mezcla in np.linspace(0, 1, 21):  # 0 = lo más cerca posible, 1 = todo parejo
        g = (1 - mezcla) * avance + mezcla * uniforme
        pos = (objetivo[0] + np.concatenate([[0], np.cumsum(g)[:-1]])) % perimetro
        idx_curva = np.searchsorted(arco, pos) % len(curva)
        Q = curva[idx_curva]
        cuerda = np.linalg.norm(Q - np.roll(Q, -1, axis=0), axis=1)
        if len(np.unique(idx_curva)) == n and cuerda.min() >= arista_min:
            break

    idx_curva = optimizar_nodos_curva(
        idx_curva, curva, [nodos_a[lazo_a], nodos_b[lazo_b]], arista_min
    )
    Q = curva[idx_curva]
    cuerda = np.linalg.norm(Q - np.roll(Q, -1, axis=0), axis=1)
    return lazo_a, lazo_b, Q, bool(cuerda.min() >= arista_min - 1e-9)


def sentido_radial(nodos, elementos, eje):
    """Signo de (normal de la cara · dirección radial respecto de un eje que pasa por
    el origen); sirve para unificar las orientaciones de las caras de un tubo."""
    eje = np.asarray(eje, dtype=float)
    p = nodos[elementos]
    normal = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
    centro = p.mean(axis=1)
    radial = centro - np.outer(centro @ eje, eje)  # le quito la componente axial
    return np.sign(np.einsum("ij,ij->i", normal, radial))


def calidad_cuadrilateros(nodos, elementos):
    """Forma de cada cuadrilátero de una malla (ver forma_cuadrilateros)."""
    p = nodos[elementos]
    return forma_cuadrilateros(p[:, 0], p[:, 1], p[:, 2], p[:, 3])


def armar_union(r_d, r_s, curva, eje_d, eje_s):
    """Une los dos tubos recortados con nodos comunes sobre la curva.

    1. Los lazos de los dos huecos se igualan en cantidad de nodos (al que tiene
       menos se le agregan nodos descartando caras del borde).
    2. Se crean los nodos comunes sobre la curva.
    3. Cada tubo se une a esos nodos con una franja de cuadriláteros.
    Devuelve un dict con la malla completa (nodos, elementos, grupo por elemento).
    """
    n_d, n_s = cantidad_de_nodos_del_lazo(r_d), cantidad_de_nodos_del_lazo(r_s)
    if n_s < n_d:
        r_s = dict(r_s, conservar=ampliar_lazo(r_s, curva, n_d))
    elif n_d < n_s:
        r_d = dict(r_d, conservar=ampliar_lazo(r_d, curva, n_s))

    nodos_d, el_d, lazo_d = lazo_final(r_d)
    nodos_s, el_s, lazo_s = lazo_final(r_s)
    arista_min = FACTOR_ARISTA_MIN * max(r_d["tam"], r_s["tam"])
    lazo_d, lazo_s, nodos_curva, cumple = nodos_comunes_en_la_curva(
        nodos_d, lazo_d, nodos_s, lazo_s, curva, arista_min
    )

    # numeración global: [grilla directo | grilla derivación | nodos de la curva]
    n_dg, n_sg, n = len(nodos_d), len(nodos_s), len(nodos_curva)
    nodos = np.vstack([nodos_d, nodos_s, nodos_curva])
    el_s = el_s + n_dg
    lazo_s = lazo_s + n_dg
    id_curva = n_dg + n_sg + np.arange(n)
    sig = (np.arange(n) + 1) % n

    def franja(lazo, el_grilla, eje):
        q = np.column_stack([lazo, lazo[sig], id_curva[sig], id_curva])
        # misma orientación (normal hacia afuera) que las caras de la grilla
        flip = sentido_radial(nodos, q, eje) != np.median(sentido_radial(nodos, el_grilla, eje))
        q[flip] = q[flip][:, ::-1]
        return q

    relleno_d = franja(lazo_d, el_d, eje_d)
    relleno_s = franja(lazo_s, el_s, eje_s)

    elementos = np.vstack([el_d, relleno_d, el_s, relleno_s])
    grupo = np.repeat([0, 1, 2, 3], [len(el_d), len(relleno_d), len(el_s), len(relleno_s)])
    q = calidad_cuadrilateros(nodos, np.vstack([relleno_d, relleno_s]))
    return dict(
        nodos=nodos,
        elementos=elementos,
        grupo=grupo,
        id_curva=id_curva,
        tam_d=r_d["tam"],
        tam_s=r_s["tam"],
        calidad=q,
        cumple=cumple,
        # las uniones que no cumplen la arista mínima solo se eligen si no hay otra
        puntaje=q.mean() + 0.5 * q.min() - (0 if cumple else 10),
    )


# =============================================================================
# 3. CURVA DE INTERSECCIÓN
# =============================================================================
curva = curva_interseccion(OD_1, OD_2, Alpha_2)
curva_local = f.Giro(curva, np.radians(-Alpha_2))  # curva en el sistema de la derivación
eje_d = (1, 0, 0)
eje_s = (np.cos(np.radians(Alpha_2)), 0, np.sin(np.radians(Alpha_2)))


# =============================================================================
# 4. TUBO DIRECTO (eje x)
# =============================================================================
nodos_d, elementos_d = f.Tube(OD_1, L_directo, N_directo)
nodos_d = f.traslacion(nodos_d, np.array([-H_1, 0, 0]))
elementos_d = elementos_d.astype(int)[:, [0, 1, 3, 2]]  # orden perimetral

# nodos del directo dentro del cilindro lateral: giro los nodos -Alpha_2 solo para
# evaluar la condición (lado de la derivación: x' > 0)
nodos_dm = f.Giro(nodos_d, np.radians(-Alpha_2))
dentro_d = (nodos_dm[:, 1] ** 2 + nodos_dm[:, 2] ** 2 < OD_2**2 / 4) & (nodos_dm[:, 0] > 0)
cruzada_d = caras_cruzadas_por_curva(
    curva, nodos_d[::N_directo, 0], N_directo, len(elementos_d)
)
recorte_d = preparar_recorte(nodos_d, elementos_d, dentro_d, cruzada_d, curva)


# =============================================================================
# 5. DERIVACIÓN (eje a Alpha_2 grados) Y UNIÓN
# =============================================================================
def recortar_derivacion(n_div):
    """Tubo de la derivación con n_div divisiones, recortado contra el directo."""
    # genero el tubo en su sistema propio (eje = x) y después lo giro a su posición
    nodos_loc, elementos = f.Tube(OD_2, L_salida, n_div)
    elementos = elementos.astype(int)[:, [0, 1, 3, 2]]  # orden perimetral
    nodos = f.Giro(nodos_loc, np.radians(Alpha_2))
    dentro = nodos[:, 1] ** 2 + nodos[:, 2] ** 2 < OD_1**2 / 4  # dentro del directo
    cruzada = caras_cruzadas_por_curva(curva_local, nodos_loc[::n_div, 0], n_div, len(elementos))
    return preparar_recorte(nodos, elementos, dentro, cruzada, curva)


# Con divisiones impares el lazo de la derivación tiene cantidad impar de nodos y no
# puede igualarse al del directo (que siempre es par), por eso solo se prueban pares.
if AJUSTAR_DIVISIONES:
    candidatas = [n for n in range(N_salida - BUSQUEDA, N_salida + BUSQUEDA + 1) if n >= 8 and n % 2 == 0]
else:
    candidatas = [N_salida]

union = None
for n_div in candidatas:
    try:
        prueba = armar_union(recorte_d, recortar_derivacion(n_div), curva, eje_d, eje_s)
    except ValueError as error:
        print(f"N_salida={n_div}: descartado ({error})")
        continue
    prueba["n_div"] = n_div
    if union is None or prueba["puntaje"] > union["puntaje"]:
        union = prueba
if union is None:
    raise RuntimeError("Ninguna cantidad de divisiones de la derivación dio una unión válida")

if not union["cumple"]:
    print(
        "AVISO: no se logró que todas las aristas sobre la curva sean >= "
        f"{FACTOR_ARISTA_MIN} x tamaño de elemento (hay demasiados nodos para su longitud)"
    )
nodos, elementos, grupo = union["nodos"], union["elementos"], union["grupo"]
N_salida = union["n_div"]  # divisiones finales de la derivación
q = union["calidad"]
print(
    f"Derivación: {N_salida} divisiones (tamaño de elemento {union['tam_s']:.1f} contra "
    f"{union['tam_d']:.1f} del directo) | {len(union['id_curva'])} nodos comunes sobre la curva\n"
    f"Relleno: {len(q)} cuadriláteros, calidad media {q.mean():.2f} y mínima {q.min():.2f} "
    f"(1 = cuadrado perfecto)"
)


# =============================================================================
# 6. GRÁFICO
# =============================================================================
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

colores_grupo = np.array(["lightsteelblue", "khaki", "palegreen", "moccasin"])
mostrar = np.ones(len(elementos), dtype=bool) if MOSTRAR_LATERAL else grupo < 2
ax.add_collection3d(
    Poly3DCollection(
        nodos[elementos[mostrar]],
        facecolors=colores_grupo[grupo[mostrar]],
        edgecolors="k",
        linewidths=0.6,
        alpha=0.9,
    )
)
usados = np.unique(elementos[mostrar])
es_curva = np.isin(usados, union["id_curva"])
ax.scatter(*nodos[usados[~es_curva]].T, color="r", s=10)  # nodos de la grilla
ax.scatter(*nodos[usados[es_curva]].T, color="g", s=25)  # nodos comunes sobre la curva
ax.plot(*np.vstack([curva, curva[:1]]).T, color="b", linewidth=1.5, label="Curva de intersección")

# límites y proporción real (x = eje del tubo directo)
mins, maxs = nodos[usados].min(axis=0), nodos[usados].max(axis=0)
ax.set_xlim(mins[0], maxs[0])
ax.set_ylim(mins[1], maxs[1])
ax.set_zlim(mins[2], maxs[2])
ax.set_box_aspect(maxs - mins)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()
