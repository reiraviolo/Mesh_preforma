"""Malla preformada de una tee por BLOQUES (solo cuadriláteros).

Idea (ver croquis): las superficies se dividen en bloques estructurados.

  * Tubo directo: anillos con N_directo nodos equiespaciados en la circunferencia
    (las dos bocas quedan equidistanciadas). La posición de los anillos a lo largo
    del eje SE ACOMODA para que haya una "ventana" rectangular alrededor de la
    derivación: x = xa ... xb y dos generatrices a cada lado.
  * Ventana: se rellena con una corona de K capas (malla en O) entre el borde de la
    ventana y la curva de intersección.
  * Banda: las columnas de la ventana son mucho más largas en x que en θ. Fuera de
    la ventana cada una se divide en DIV_BANDA; las dos filas que tocan el borde
    llevan una plantilla en ROMBO, que entra con 1 arista del lado de la ventana y
    sale con DIV_BANDA del otro, de modo que el borde no gana nodos y la derivación
    conserva sus divisiones.
  * Corona: el número de capas varía según el grosor local; los cambios (de a 2) se
    resuelven con la misma plantilla en ROMBO, girada 90°.
  * Derivación: N_salida nodos equiespaciados en la circunferencia (boca equidistanciada).
    Sus generatrices arrancan EXACTAMENTE en los nodos de la curva, que son los mismos
    nodos de la corona del directo (coincidentes por construcción).

Pasos:  1. Parámetros  2. Curva de intersección  3. Tubo directo (grilla + ventana)
        3.5 Banda refinada (rombos)  4. Corona  5. Derivación  6. Control  7. Gráfico
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from collections import Counter

# =============================================================================
# 1. PARÁMETROS
# =============================================================================
OD_1 = 200  # OD tubo directo
OD_2 = 200  # OD derivación (condición: OD_2 <= OD_1)
H_1 = 200  # de entrada a intersección
H_2 = 200  # de intersección a salida lateral
H_3 = 300  # de intersección a salida final
Alpha_2 = 60  # ángulo de la derivación [grados]

N_directo = 100  # divisiones de la boca del directo (se redondea a múltiplo de 4)
CAPAS = 1  # capas mínimas de elementos entre la ventana y la curva (corona)
MARGEN = 2  # separación mínima entre la curva y el borde de la ventana, en tamaños de elemento
CAPAS_TOPES = (None, 7, 5, 3)  # topes de capas de la corona que prueba la autocalibración
                               # (None = sin tope); se arma con cada uno y se queda con el mejor
JAC_MIN = 0.2  # jacobiano escalado mínimo aceptable (umbral habitual en elementos finitos)
SUAVIZADO = 100  # iteraciones de suavizado: la grilla se deforma para adaptarse a la curva
DIST_MIN = 1.0  # los nodos del borde de la ventana no se acercan a la curva a menos de DIST_MIN·tamaño
BANDA = 8  # anillos a cada lado de la ventana que pueden deformarse
TRANSICION = 3.0  # largo de la zona de transición de tamaños junto a la ventana (en tamaños de ventana)
N_X = None  # celdas de la ventana en x (None = las mínimas: N_directo/10, así hay pocos nodos agregados)
DIV_BANDA = None  # divisiones en x de las columnas, fuera de la ventana (None = automatico; solo 1 o 3)
LONG_MIN = 0.4  # ningún nodo se acerca a un vecino a menos de LONG_MIN·tamaño (evita celdas colapsadas en las puntas)

R1, R2 = OD_1 / 2, OD_2 / 2
alfa = np.radians(Alpha_2)
eje_s = np.array([np.cos(alfa), 0.0, np.sin(alfa)])  # eje de la derivación
perp_s = np.array([-np.sin(alfa), 0.0, np.cos(alfa)])  # transversal (plano xz)

N_directo = max(8, 4 * int(round(N_directo / 4)))
dth = 2 * np.pi / N_directo
tam = R1 * dth  # tamaño de elemento del directo


# =============================================================================
# 2. CURVA DE INTERSECCIÓN  (φ = ángulo en la sección de la derivación)
# =============================================================================
def punto_curva(phi):
    """Punto de la curva para el ángulo φ de la derivación: y=R2 cosφ, v=R2 sinφ."""
    phi = np.asarray(phi, dtype=float)
    y = R2 * np.cos(phi)
    z = np.sqrt(R1**2 - y**2)
    x = (z * np.cos(alfa) - R2 * np.sin(phi)) / np.sin(alfa)
    return np.stack([x, y, z], axis=-1)


fi_denso = np.linspace(0, 2 * np.pi, 4000, endpoint=False)
p_denso = punto_curva(fi_denso)
x_min, x_max = p_denso[:, 0].min(), p_denso[:, 0].max()
th_min = np.arctan2(p_denso[:, 2], p_denso[:, 1]).min()
th_max = np.arctan2(p_denso[:, 2], p_denso[:, 1]).max()

# =============================================================================
# 3. TUBO DIRECTO: grilla de anillos con una ventana rectangular
# =============================================================================
# 3.0 N_directo mínimo para que la ventana entre en el tubo (si hace falta se sube de a 4)
def ventana(N):
    dth_, tam_ = 2 * np.pi / N, R1 * 2 * np.pi / N
    n_th_ = 2 * int(np.ceil((np.pi / 2 - th_min) / dth_ + CAPAS))
    ancho_ = x_max - x_min + 2 * max(CAPAS, MARGEN) * tam_
    n_x_nat = int(np.ceil(ancho_ / tam_))
    n_x_ = n_x_nat + n_x_nat % 2  # paso parejo: las celdas de la ventana miden como las del directo
    return n_th_, max(4, n_x_)


N_pedido = N_directo
while ventana(N_directo)[0] > N_directo - 2:
    N_directo += 4
if N_directo != N_pedido:
    print(f"AVISO: N_directo={N_pedido} es muy chico para esta derivación; se usa N_directo={N_directo}")
dth = 2 * np.pi / N_directo
tam = R1 * dth

# 3.1 columnas de la ventana (simétricas respecto de la cúspide θ = 90°)
n_th = 2 * int(np.ceil((np.pi / 2 - th_min) / dth + CAPAS))  # celdas en θ (par)
if n_th > N_directo - 2:
    raise ValueError("N_directo demasiado chico para la derivación")
j1 = N_directo // 4 - n_th // 2
j2 = N_directo // 4 + n_th // 2

# 3.2 anillos. La ventana tiene n_x celdas (par) y cubre la curva + margen, con el mismo
#     paso que el resto del directo. Las divisiones de la derivación no salen de acá:
#     la corona las ajusta con rombos (sección 4).
margen = max(CAPAS, MARGEN) * tam  # separación entre la curva y el borde de la ventana
n_x_natural = int(np.ceil((x_max - x_min + 2 * margen) / tam))
n_x = ventana(N_directo)[1] if N_X is None else N_X
if n_x % 2 or n_x < 4:
    raise ValueError("N_X debe ser par y >= 4")
x_c = (x_min + x_max) / 2
ancho = x_max - x_min + 2 * margen
tam_x = ancho / n_x  # tamaño de las celdas de la ventana en x (>= tam si se limitó)
xa, xb = x_c - ancho / 2, x_c + ancho / 2


def anillos(a, b, s_a, s_b):
    """Posiciones entre a y b con tamaño de celda que varía linealmente de s_a a s_b."""
    L = b - a
    if abs(s_b - s_a) < 1e-9 * s_a:
        n_real = L / s_a
        k = 0.0
    else:
        k = (s_b - s_a) / L
        n_real = np.log(s_b / s_a) / k
    n = max(1, int(round(n_real)))
    u = np.linspace(0, n_real, n + 1)
    pos = u / n_real if k == 0 else (np.exp(k * u) - 1) * s_a / (s_b - s_a)
    return a + L * pos / pos[-1]


def anillos_con_rampa(a, b, s_a, s_b, rampa_al_final):
    """Anillos entre a y b: tamaño constante (el de la boca) y una rampa corta junto a la ventana."""
    largo = b - a
    L = min(TRANSICION * max(s_a, s_b), 0.8 * largo)
    if rampa_al_final:  # la ventana está en b
        r = anillos(b - L, b, tam, tam_x)
        m = max(1, int(round((largo - L) / tam)))
        return np.concatenate([np.linspace(a, b - L, m + 1)[:-1], r])
    r = anillos(a, a + L, tam_x, tam)  # la ventana está en a
    m = max(1, int(round((largo - L) / tam)))
    return np.concatenate([r[:-1], np.linspace(a + L, b, m + 1)])


if DIV_BANDA is None:  # se divide en 3 solo si asi las celdas quedan mas cuadradas
    DIV_BANDA = 3 if tam_x / tam > np.sqrt(3) else 1
tam_x_ef = tam_x / DIV_BANDA  # tamano real en x fuera de la ventana, ya dividido por los rombos
X_izq = anillos_con_rampa(-H_1, xa, tam, tam_x_ef, True)
X_der = anillos_con_rampa(xb, H_3, tam_x_ef, tam, False)
X = np.concatenate([X_izq[:-1], np.linspace(xa, xb, n_x + 1)[:-1], X_der])
n_izq = len(X_izq) - 1
i1, i2 = n_izq, n_izq + n_x  # anillos de la ventana
n_anillos = len(X)


def id_grilla(i, j):
    return i * N_directo + (j % N_directo)


theta = np.arange(N_directo) * dth
nodos_grilla = np.array(
    [[x, R1 * np.cos(t), R1 * np.sin(t)] for x in X for t in theta]
)
cuads_grilla = []
for i in range(n_anillos - 1):
    for j in range(N_directo):
        dentro = i1 <= i < i2 and (j - j1) % N_directo < n_th
        banda = i1 <= i < i2 and not dentro  # se arma aparte, con rombos
        if not dentro and not (banda and DIV_BANDA > 1):
            cuads_grilla.append(
                [id_grilla(i, j), id_grilla(i, j + 1), id_grilla(i + 1, j + 1), id_grilla(i + 1, j)]
            )
cuads_grilla = np.array(cuads_grilla)

# =============================================================================
# 3.5 BANDA REFINADA (rombos)
# =============================================================================
# Las columnas de la ventana son tam_x (mucho mas largas en x que en theta). Fuera
# de la ventana cada una se divide en DIV_BANDA, pero el borde de la ventana no puede
# ganar nodos (fijan N_s). Por eso las dos filas que tocan el borde llevan una
# plantilla en rombo: entran con 1 arista del lado de la ventana y salen con
# DIV_BANDA del otro lado, agregando 2 elementos. Topologicamente hace falta que
# DIV_BANDA sea impar (1 + 1 + DIV_BANDA + 1 debe ser par).
n_estruct = len(nodos_grilla)
if DIV_BANDA not in (1, 3):
    raise ValueError("DIV_BANDA debe ser 1 o 3")
nodos_extra, cuads_banda = [], []
if DIV_BANDA > 1:
    def _nodo(x, th):
        nodos_extra.append([x, R1 * np.cos(th), R1 * np.sin(th)])
        return n_estruct + len(nodos_extra) - 1

    n_banda = N_directo - n_th  # filas de celdas de la banda
    filas = [(j2 + r) % N_directo for r in range(n_banda + 1)]  # de un borde al otro
    for i in range(i1, i2):
        xs = [X[i] + (X[i + 1] - X[i]) * k / DIV_BANDA for k in range(1, DIV_BANDA)]
        # nodos de las lineas interiores (las dos filas del borde no llevan)
        M = {r: [_nodo(x, filas[r] * dth) for x in xs] for r in range(1, n_banda)}
        for r in range(n_banda):
            A, B = id_grilla(i, filas[r]), id_grilla(i + 1, filas[r])
            D, C = id_grilla(i, filas[r + 1]), id_grilla(i + 1, filas[r + 1])
            if r == 0:  # rombo: 1 arista arriba (borde) -> DIV_BANDA abajo
                th_p = (filas[0] + 0.5) * dth
                P = [_nodo(x, th_p) for x in xs]
                cuads_banda += [[A, B, P[1], P[0]], [A, P[0], M[1][0], D],
                                [P[0], P[1], M[1][1], M[1][0]], [P[1], B, C, M[1][1]]]
            elif r == n_banda - 1:  # rombo: DIV_BANDA arriba -> 1 arista abajo (borde)
                th_q = (filas[r] + 0.5) * dth
                Q = [_nodo(x, th_q) for x in xs]
                t = M[r]
                cuads_banda += [[Q[0], Q[1], C, D], [A, t[0], Q[0], D],
                                [t[0], t[1], Q[1], Q[0]], [t[1], B, C, Q[1]]]
            else:  # fila regular: DIV_BANDA celdas
                t, b = M[r], M[r + 1]
                cuads_banda += [[A, t[0], b[0], D], [t[0], t[1], b[1], b[0]], [t[1], B, C, b[1]]]
nodos_grilla = np.vstack([nodos_grilla, np.array(nodos_extra).reshape(-1, 3)])
cuads_grilla = np.array(cuads_grilla.tolist() + cuads_banda)

# =============================================================================
# ANILLO DE CAPAS VARIABLES (rombos)
# =============================================================================
# Malla entre dos lazos cerrados de N nodos cuando el número de capas tiene que
# cambiar a lo largo del lazo (porque el ancho a cubrir no es parejo). Donde cambia
# se usa la plantilla en ROMBO: entra 1 arista por un lado y salen 3 por el otro,
# agregando 2 elementos. Por eso el número de capas solo puede cambiar de a 2: si
# cambiara de a 1 el parche tendría un número impar de aristas de borde, imposible
# con cuadriláteros. Se usa para la corona y para la derivación.
def capas_por_tamano(ancho, paso, c_min, c_max=None):
    """Capas en cada posición para que el tamaño a lo ancho se parezca a `paso`.
    Todas de la misma paridad y con saltos de 0 o 2 entre posiciones vecinas."""
    tope = len(ancho) if c_max is None else max(0, (c_max - c_min) // 2)
    capas = c_min + 2 * np.clip(np.round((ancho / paso - c_min) / 2).astype(int), 0, tope)
    for _ in range(len(capas)):
        recortado = np.minimum(capas, np.minimum(np.roll(capas, 1), np.roll(capas, -1)) + 2)
        if np.array_equal(recortado, capas):
            break
        capas = recortado
    return capas


def nodos_de_anillo(capas):
    """Nodos que agrega `anillo_variable` (sin contar los dos lazos)."""
    return int((capas - 1).sum()) + 2 * int((capas != np.roll(capas, -1)).sum())


def anillo_variable(ids_ext, ids_int, capas, punto, proyecta, base):
    """Cuadriláteros entre el lazo exterior y el interior con `capas[k]` capas en cada
    posición. `punto(k, v)` da la posición en k a la fracción v (0 = exterior, 1 = interior)
    y `proyecta(p)` devuelve p sobre la superficie. `base` es el primer id libre.
    Devuelve (coordenadas de los nodos nuevos, cuadriláteros)."""
    N = len(ids_ext)
    ini = base + np.concatenate([[0], np.cumsum(capas - 1)])[:-1]
    base_rombo = base + int((capas - 1).sum())
    coord = np.zeros((nodos_de_anillo(capas), 3))
    fila = []
    for k in range(N):
        fila.append([ids_ext[k]] + [ini[k] + j for j in range(capas[k] - 1)] + [ids_int[k]])
        for j in range(1, capas[k]):
            coord[ini[k] - base + j - 1] = punto(k, j / capas[k])
    cuads, usados = [], 0
    for k in range(N):
        kk = (k + 1) % N
        L, R, c0, c1 = fila[k], fila[kk], int(capas[k]), int(capas[kk])
        if c0 == c1:
            cuads += [[L[j], R[j], R[j + 1], L[j + 1]] for j in range(c0)]
            continue
        if abs(c1 - c0) != 2:
            raise RuntimeError("el número de capas solo puede cambiar de a 2")
        m = min(c0, c1) // 2  # capa donde nacen (o mueren) las dos capas nuevas
        P = [base_rombo + usados, base_rombo + usados + 1]
        usados += 2
        if c1 == c0 + 2:  # 1 arista en la posición k -> 3 en la posición k+1
            A, B, M1, M2, C, D = L[m], R[m], R[m + 1], R[m + 2], R[m + 3], L[m + 1]
            for j in (0, 1):
                coord[P[j] - base] = proyecta(
                    0.55 * punto(k, (m + (j + 1) / 3) / c0) + 0.45 * punto(kk, (m + 1 + j) / c1)
                )
            cuads += [[L[j], R[j], R[j + 1], L[j + 1]] for j in range(m)]
            cuads += [[L[j], R[j + 2], R[j + 3], L[j + 1]] for j in range(m + 1, c0)]
        else:  # 3 aristas en la posición k -> 1 en la posición k+1
            A, B, M1, M2, C, D = R[m], L[m], L[m + 1], L[m + 2], L[m + 3], R[m + 1]
            for j in (0, 1):
                coord[P[j] - base] = proyecta(
                    0.55 * punto(kk, (m + (j + 1) / 3) / c1) + 0.45 * punto(k, (m + 1 + j) / c0)
                )
            cuads += [[L[j], R[j], R[j + 1], L[j + 1]] for j in range(m)]
            cuads += [[L[j], R[j - 2], R[j - 1], L[j + 1]] for j in range(m + 3, c0)]
        cuads += [[A, B, M1, P[0]], [P[0], M1, M2, P[1]], [P[1], M2, C, D], [A, P[0], P[1], D]]
    return coord, np.array(cuads)


def sobre_directo(p):
    """Devuelve p sobre el cilindro del directo."""
    r = np.hypot(p[1], p[2])
    return np.array([p[0], p[1] * R1 / r, p[2] * R1 / r])


def sobre_derivacion(p):
    """Devuelve p sobre el cilindro de la derivación."""
    t = p @ eje_s
    rv = p - t * eje_s
    return t * eje_s + rv * R2 / np.linalg.norm(rv)


# =============================================================================
# 4. CORONA (malla en O) entre la ventana y la curva
# =============================================================================
# Borde de la ventana recorrido en el mismo sentido en que φ crece sobre la curva:
#   lado "abajo" (θ=θ1, x decrece) -> lado izquierdo (θ crece) -> lado "arriba" (x crece)
#   -> lado derecho (θ decrece).  Empieza en la mitad de n_x/2 antes de φ = 0.
# abajo: de (i2,j1) hasta (i1,j1), pero el índice q=0 está en la mitad
abajo = [(i2 - k, j1) for k in range(n_x)]
izq = [(i1, j1 + k) for k in range(n_th)]
arriba = [(i1 + k, j2) for k in range(n_x)]
der = [(i2, j2 - k) for k in range(n_th)]
lazo = abajo + izq + arriba + der  # sentido: φ creciente, empezando en q = -n_x/2
N_borde = len(lazo)  # nodos del borde de la ventana
lazo = lazo[n_x // 2 :] + lazo[: n_x // 2]  # ahora el nodo 0 corresponde a φ = 0
assert N_borde == 2 * (n_x + n_th)

# Ángulos φ de las 4 esquinas de la ventana: se unen con los puntos de la curva que
# están en la misma dirección desde el centro de la ventana.
pl_denso = np.stack([p_denso[:, 0], R1 * np.arctan2(p_denso[:, 2], p_denso[:, 1])], axis=1)
centro_v = np.array([x_c, R1 * np.pi / 2])
esquinas = [(xa, j1), (xa, j2), (xb, j2), (xb, j1)]  # sentido de φ creciente
phi_esq = []
for xe, je in esquinas:
    dirv = np.array([xe, R1 * je * dth]) - centro_v
    ang = np.arctan2(*(pl_denso - centro_v)[:, ::-1].T) - np.arctan2(dirv[1], dirv[0])
    phi_esq.append(fi_denso[np.argmin(np.abs(np.angle(np.exp(1j * ang))))])
phi_esq = np.unwrap(np.array(phi_esq + [phi_esq[0] + 2 * np.pi]))
phi_esq = np.array([phi_esq[0], phi_esq[1], np.pi, phi_esq[2], phi_esq[3], 2 * np.pi, phi_esq[4]])
phi_ext = np.concatenate([fi_denso, fi_denso + 2 * np.pi])
pl_ext = np.vstack([pl_denso, pl_denso])
arco = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pl_ext, axis=0), axis=1))])
arco_esq = np.array([np.interp(f, phi_ext, arco) for f in phi_esq])

# ---- cuántos nodos lleva la curva ------------------------------------------
# El borde de la ventana tiene el paso del directo. La curva es bastante más corta que
# el borde, así que para el mismo tamaño de elemento necesita menos nodos; si se usaran
# los del borde, la derivación entera saldría demasiado fina. La diferencia la absorbe
# un anillo de rombos en la corona (3 aristas del borde entran en 1 de la curva).
# El lazo se parte en 6 segmentos entre esquinas: θ, x/2, x/2, θ, x/2, x/2.
nb_seg = np.array([n_th, n_x // 2, n_x // 2, n_th, n_x // 2, n_x // 2])
arco_seg = np.diff(arco_esq)
a_th = (arco_seg[0] + arco_seg[3]) / 2  # arco de un lado en θ (por simetría, los 2 iguales)
a_x = (arco_seg[1] + arco_seg[2] + arco_seg[4] + arco_seg[5]) / 4  # arco de medio lado recto


def nodos_de_segmento(ideal, nb):
    """Aristas de la curva en un segmento con nb del borde: entre ceil(nb/3) (cada celda
    toma a lo sumo 3) y nb, y con la paridad que hace par a la diferencia."""
    piso = -(-nb // 3)
    piso += (nb - piso) % 2
    c = int(round(ideal))
    c -= (nb - c) % 2
    return int(min(max(c, piso), nb))


n_obj = (2 * a_th + 4 * a_x) / tam  # nodos que pide el tamaño de elemento objetivo
c_th = nodos_de_segmento(n_obj * a_th / (2 * a_th + 4 * a_x), n_th)
c_x = nodos_de_segmento(n_obj * a_x / (2 * a_th + 4 * a_x), n_x // 2)
nc_seg = np.array([c_th, c_x, c_x, c_th, c_x, c_x])
N_s = int(nc_seg.sum())  # divisiones de la curva y de la derivación

# φ de los nodos de la curva: entre esquinas se reparten por longitud de arco (en el
# desarrollo x, R·θ), así los elementos de la corona quedan lo menos deformados posible.
pc_esq = np.concatenate([[c_x], c_x + np.cumsum(nc_seg)])
phi_c = np.zeros(N_s)
for k in range(6):
    for p in range(pc_esq[k], pc_esq[k + 1]):
        t = (p - pc_esq[k]) / (pc_esq[k + 1] - pc_esq[k])
        phi_c[p % N_s] = np.interp(arco_esq[k] + (arco_esq[k + 1] - arco_esq[k]) * t, arco, phi_ext) % (2 * np.pi)
nodos_curva = punto_curva(phi_c)

# ángulos φ de la boca de la derivación: equiespaciados, alineados con los de la curva
phi_boca = np.arange(N_s) * 2 * np.pi / N_s
desfase = np.angle(np.mean(np.exp(1j * (phi_c - phi_boca))))
phi_boca = phi_boca + desfase

id_grilla_borde = np.array([id_grilla(i, j) for i, j in lazo])
n_g = len(nodos_grilla)


def a_plano(p):
    """(x, θ) de puntos sobre el cilindro directo."""
    return np.stack([p[..., 0], np.arctan2(p[..., 2], p[..., 1])], axis=-1)


def a_3d(pl):
    """(x, θ) -> punto sobre el cilindro directo."""
    return np.stack([pl[..., 0], R1 * np.cos(pl[..., 1]), R1 * np.sin(pl[..., 1])], axis=-1)


a_borde = np.array([[X[i], j * dth] for i, j in lazo])  # (x, θ) con θ continuo (j puede salirse de 0..N-1)
a_curva = a_plano(nodos_curva)

# ---- anillo de transición: N_borde aristas -> N_s ---------------------------
# Cada celda toma 1 arista del borde (celda normal) o 3 (rombo). Los rombos se reparten
# parejos dentro de cada segmento, tantos como (nb - nc) / 2.
# Un rombo toma 3 aristas del borde, así que necesita lugar a lo ancho: se ponen donde
# la corona es más gruesa (nunca dos pegados), no reparados a ciegas.
pb_esq = np.concatenate([[n_x // 2], n_x // 2 + np.cumsum(nb_seg)])
prov = np.zeros(N_s, dtype=int)  # mapeo provisorio nodo de la curva -> nodo del borde
for k in range(6):
    for j in range(nc_seg[k]):
        prov[(pc_esq[k] + j) % N_s] = (pb_esq[k] + round(j * nb_seg[k] / nc_seg[k])) % N_borde
grosor_ref = np.linalg.norm(  # primera estimación; la autocalibración la reemplaza por la real
    np.stack([a_borde[prov, 0], R1 * a_borde[prov, 1]], axis=1)
    - np.stack([a_curva[:, 0], R1 * a_curva[:, 1]], axis=1), axis=1
)

def suavizar(nodos, cuads, libres_ids, iteraciones, w=0.6, cuidar=(), curva=None, dist_min=0.0, long_min=0.0,
             proyecta=None):
    """Laplaciano sobre la superficie; los nodos de `libres_ids` se mueven, el resto no."""
    proyecta = proyecta or sobre_directo
    vecinos = [set() for _ in range(len(nodos))]
    for q in cuads:
        for a in range(4):
            vecinos[q[a]].update((q[(a + 1) % 4], q[(a - 1) % 4]))
    nodos = nodos.copy()
    for _ in range(iteraciones):
        nuevo = nodos.copy()
        for i in libres_ids:
            if not vecinos[i]:
                continue
            m = nodos[list(vecinos[i])].mean(axis=0)
            q = proyecta(nodos[i] + w * (m - nodos[i]))  # de nuevo sobre la superficie
            if i in cuidar and np.linalg.norm(curva - q, axis=1).min() < dist_min:
                continue  # un nodo del borde de la ventana no se acerca demasiado a la curva
            if long_min and np.linalg.norm(nodos[list(vecinos[i])] - q, axis=1).min() < long_min:
                continue  # no se colapsan aristas: el nodo no se acerca demasiado a ningún vecino
            nuevo[i] = q
        nodos = nuevo
    return nodos




def calidad(nodos, elementos):
    """(dispersión de tamaño p95/p5, jacobiano escalado mínimo) de una malla."""
    p = nodos[elementos]
    lado = np.sqrt(0.5 * np.linalg.norm(np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1]), axis=1))
    jac = np.min([
        np.linalg.norm(np.cross(p[:, (a + 1) % 4] - p[:, a], p[:, (a - 1) % 4] - p[:, a]), axis=1)
        / np.maximum(np.linalg.norm(p[:, (a + 1) % 4] - p[:, a], axis=1)
                     * np.linalg.norm(p[:, (a - 1) % 4] - p[:, a], axis=1), 1e-12)
        for a in range(4)
    ])
    return np.percentile(lado, 95) / np.percentile(lado, 5), jac


# AUTOCALIBRACIÓN. El ancho de la corona antes de suavizar no sirve para repartir las
# capas: el suavizado mueve bastante el borde de la ventana. Y cuántas capas conviene
# poner depende de la geometría, así que no hay un tope bueno para todos los casos. Se
# arma la corona con cada tope de CAPAS_TOPES, se suaviza, se mide, y se guarda la mejor:
# primero que el jacobiano sea usable, después la uniformidad de tamaño. Desde la segunda
# pasada el reparto usa el ancho real medido en la anterior.
mejor = None
for _pasada, CAPAS_MAX in enumerate(CAPAS_TOPES):
    tipo = []
    for k in range(6):
        mm, nc = (nb_seg[k] - nc_seg[k]) // 2, nc_seg[k]
        t = [1] * nc
        # peso de cada celda: el grosor disponible en sus dos extremos
        peso = [min(grosor_ref[(pc_esq[k] + j) % N_s], grosor_ref[(pc_esq[k] + j + 1) % N_s]) for j in range(nc)]
        elegidas = []
        for j in sorted(range(nc), key=lambda q: -peso[q]):
            if len(elegidas) == mm:
                break
            if all(abs(j - e) > 1 for e in elegidas):
                elegidas.append(j)
        for j in sorted(range(nc), key=lambda q: -peso[q]):  # si no alcanzó, se relaja la separación
            if len(elegidas) == mm:
                break
            if j not in elegidas:
                elegidas.append(j)
        for j in elegidas:
            t[j] = 3
        tipo += t
    ini_cel, s_b = [], n_x // 2  # el segmento 0 empieza en el nodo n_x/2 del borde
    for q in range(N_s):
        ini_cel.append(s_b % N_borde)
        s_b += tipo[q]
    assert s_b - n_x // 2 == N_borde
    borde_de = np.zeros(N_s, dtype=int)  # nodo del borde del que cuelga cada nodo del anillo
    for q in range(N_s):
        borde_de[(c_x + q) % N_s] = ini_cel[q]

    # ---- capas variables entre el anillo y la curva -----------------------------
    # El grosor de la corona cambia mucho alrededor del lazo (mínimo en el medio de cada
    # lado, máximo en las esquinas de la ventana), así que con un número fijo de capas los
    # elementos de las esquinas quedan mucho más gruesos que largos.
    pl_b = np.stack([a_borde[borde_de, 0], R1 * a_borde[borde_de, 1]], axis=1)
    pl_c = np.stack([a_curva[:, 0], R1 * a_curva[:, 1]], axis=1)
    paso = (np.linalg.norm(np.roll(pl_b, -1, 0) - pl_b, axis=1) + np.linalg.norm(np.roll(pl_c, -1, 0) - pl_c, axis=1)) / 2
    capas = capas_por_tamano(grosor_ref, (paso + np.roll(paso, 1)) / 2, CAPAS, CAPAS_MAX)
    prof = 1.0 / (capas + 1)  # el anillo de transición se lleva una capa
    a_anillo = (1 - prof)[:, None] * a_borde[borde_de] + prof[:, None] * a_curva

    id_anillo = n_g + np.arange(N_s)
    base_rombo_anillo = n_g + N_s
    n_rombo_anillo = 2 * sum(1 for t in tipo if t == 3)
    base_capas = base_rombo_anillo + n_rombo_anillo
    n_corona_int = N_s + n_rombo_anillo + nodos_de_anillo(capas)
    id_curva = n_g + n_corona_int + np.arange(N_s)  # numeración de los nodos de la curva

    cuads_anillo, coord_rombo = [], []
    for q in range(N_s):
        k, k2 = (c_x + q) % N_s, (c_x + q + 1) % N_s
        o = [id_grilla_borde[(ini_cel[q] + t) % N_borde] for t in range(tipo[q] + 1)]
        if tipo[q] == 1:
            cuads_anillo.append([o[0], o[1], id_anillo[k2], id_anillo[k]])
            continue
        P = [base_rombo_anillo + len(coord_rombo), base_rombo_anillo + len(coord_rombo) + 1]
        medio = (a_anillo[k] + a_anillo[k2]) / 2
        for j in (1, 2):
            coord_rombo.append(a_3d(0.6 * a_borde[(ini_cel[q] + j) % N_borde] + 0.4 * medio))
        A, B, M1, M2, C, D = id_anillo[k], o[0], o[1], o[2], o[3], id_anillo[k2]
        cuads_anillo += [[A, B, M1, P[0]], [P[0], M1, M2, P[1]], [P[1], M2, C, D], [A, P[0], P[1], D]]

    coord_capas, cuads_capas = anillo_variable(
        id_anillo, id_curva, capas,
        lambda k, v: a_3d((1 - v) * a_anillo[k] + v * a_curva[k]), sobre_directo, base_capas,
    )
    nodos_int = np.vstack([a_3d(a_anillo), np.array(coord_rombo).reshape(-1, 3), coord_capas])
    cuads_corona = np.vstack([np.array(cuads_anillo), cuads_capas])

    # =============================================================================
    # 5. DERIVACIÓN: generatrices desde los nodos de la curva hasta la boca (t = H_2)
    # =============================================================================
    t_curva = nodos_curva @ eje_s  # posición axial de cada nodo de la curva
    tam_s = 2 * np.pi * R2 / N_s  # tamaño circunferencial de la derivación
    tam_axial = tam_s  # elementos cuadrados en la derivación (mismo tamaño axial que circunferencial)
    dphi = np.angle(np.exp(1j * (phi_boca - phi_c)))  # giro de las generatrices hasta la boca
    # Las generatrices tienen largos muy distintos (la curva sube y baja), así que con un
    # número fijo de capas los elementos de las generatrices largas quedan mucho más
    # estirados. Cada una lleva las capas que le corresponden, con rombos en los cambios.
    largo_gen = H_2 - t_curva
    if largo_gen.min() <= 0:
        raise ValueError(f"H_2={H_2:g} es demasiado corto: la boca corta la curva; usar H_2 > {t_curva.max():.0f}")
    M_min = max(1, int(round(largo_gen.min() / tam_axial)))
    capas_der = capas_por_tamano(largo_gen, tam_axial, M_min, M_min + 2 * len(largo_gen))
    M = int(np.median(capas_der))


    def punto_der(k, v):
        """Punto de la generatriz k a la fracción v (0 = curva, 1 = boca)."""
        phi = phi_c[k] + dphi[k] * v  # la última capa (boca) queda equiespaciada y plana en t = H_2
        t = t_curva[k] + largo_gen[k] * v
        return t * eje_s + R2 * np.cos(phi) * np.array([0, 1.0, 0]) + R2 * np.sin(phi) * perp_s


    id_boca = n_g + n_corona_int + N_s + np.arange(N_s)  # nodos de la boca de la derivación
    nodos_boca = np.array([punto_der(k, 1.0) for k in range(N_s)])
    nodos_der, cuads_der = anillo_variable(
        id_curva, id_boca, capas_der, punto_der, sobre_derivacion, n_g + n_corona_int + 2 * N_s,
    )

    # ---- ensamble --------------------------------------------------------------
    nodos = np.vstack([nodos_grilla, nodos_int, nodos_curva, nodos_boca, nodos_der])
    elementos = np.vstack([cuads_grilla, cuads_corona, cuads_der])
    grupo = np.concatenate(
        [np.zeros(len(cuads_grilla), int), np.ones(len(cuads_corona), int), np.full(len(cuads_der), 2)]
    )


    # ---- suavizado: la grilla del directo se deforma para adaptarse a la curva -------------
    cuads_directo = np.vstack([cuads_grilla, cuads_corona])
    anillo_de = np.arange(n_estruct) // N_directo
    libres_ids = [
        i for i in range(n_estruct) if i1 - BANDA <= anillo_de[i] <= i2 + BANDA and 0 < anillo_de[i] < n_anillos - 1
    ]
    libres_ids += list(range(n_estruct, n_g))  # nodos agregados por los rombos
    borde_ventana = set(id_grilla_borde.tolist())
    libres_ids = libres_ids + list(range(n_g, n_g + n_corona_int))  # capas intermedias de la corona
    if SUAVIZADO:
        nodos = suavizar(
            nodos, cuads_directo, libres_ids, SUAVIZADO, cuidar=borde_ventana, curva=nodos_curva, dist_min=DIST_MIN * tam,
            long_min=LONG_MIN * tam,
        )

    # los nodos interiores de la derivación también se acomodan (la curva y la boca no se mueven)
    if SUAVIZADO:
        nodos = suavizar(
            nodos, cuads_der, list(range(n_g + n_corona_int + 2 * N_s, len(nodos))), SUAVIZADO,
            proyecta=sobre_derivacion,
        )

    # ---- autocalibración --------------------------------------------------
    # El grosor de antes de suavizar no sirve para repartir las capas: el suavizado
    # mueve bastante el borde de la ventana. Con la malla ya suavizada se mide el
    # ancho real de la corona y se vuelve a repartir; en general dos pasadas alcanzan.
    ancho_real = np.linalg.norm(nodos[id_anillo] - nodos[id_curva], axis=1) / (1 - prof)
    # se guarda la mejor pasada: primero que el jacobiano sea usable, después la uniformidad
    _disp, _jac = calidad(nodos, elementos)
    _puntaje = (_jac >= JAC_MIN, -round(_disp, 2), _jac)  # en empate de uniformidad gana la calidad
    if mejor is None or _puntaje > mejor[0]:
        mejor = (_puntaje, nodos, elementos, grupo, id_curva, capas)
    grosor_ref = ancho_real  # las pasadas siguientes usan el ancho real, no el estimado



_, nodos, elementos, grupo, id_curva, capas = mejor  # se queda con la mejor pasada

# ---- orientación: normales hacia afuera de cada tubo ------------------------
def orientar(nodos, elementos, grupo):
    elementos = elementos.copy()
    p = nodos[elementos]
    normal = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
    centro = p.mean(axis=1)
    radial = centro.copy()
    radial[:, 0] = 0  # directo: radial respecto del eje x
    rs = centro - (centro @ eje_s)[:, None] * eje_s  # derivación: respecto de su eje
    radial[grupo == 2] = rs[grupo == 2]
    flip = np.einsum("ij,ij->i", normal, radial) < 0
    elementos[flip] = elementos[flip][:, ::-1]
    return elementos, flip


elementos, invertidos = orientar(nodos, elementos, grupo)

# ---- se descartan los nodos de la grilla dentro de la ventana (no pertenecen a ningún elemento)
usados = np.unique(elementos)
nueva_id = -np.ones(len(nodos), dtype=int)
nueva_id[usados] = np.arange(len(usados))
nodos = nodos[usados]
elementos = nueva_id[elementos]
id_curva = nueva_id[id_curva]

# =============================================================================
# 6. CONTROL
# =============================================================================
p = nodos[elementos]
lados = np.linalg.norm(p - np.roll(p, -1, axis=1), axis=2)
area = 0.5 * np.linalg.norm(np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1]), axis=1)
forma = 4 * area / (lados**2).sum(axis=1)  # 1 = cuadrado perfecto
aristas_ = np.sort(np.vstack([elementos[:, [a, (a + 1) % 4]] for a in range(4)]), axis=1)
cuenta = Counter(map(tuple, aristas_))
libres = sum(v == 1 for v in cuenta.values())
print(
    f"Directo: {N_directo} divisiones, {n_anillos} anillos (tamaño {tam:.1f}, en la ventana {tam_x:.1f} en x) | "
    f"Derivación: {N_s} divisiones, {M} capas (tamaño {tam_s:.1f})\n"
    f"Nodos: {len(nodos)} | Elementos: {len(elementos)} | "
    f"forma media {forma.mean():.2f}, mínima {forma.min():.2f} (1 = cuadrado)\n"
    f"Aristas en >2 caras: {sum(v > 2 for v in cuenta.values())} | "
    f"libres: {libres} (esperado {2 * N_directo + N_s})"
)

_lado = np.sqrt(
    0.5 * np.linalg.norm(np.cross(nodos[elementos][:, 2] - nodos[elementos][:, 0],
                                  nodos[elementos][:, 3] - nodos[elementos][:, 1]), axis=1)
)
print(
    f"Tamaño de elemento: objetivo {tam:.1f}, mediana {np.median(_lado):.1f}, "
    f"p5 {np.percentile(_lado, 5):.1f}, p95 {np.percentile(_lado, 95):.1f} "
    f"(dispersión p95/p5 {np.percentile(_lado, 95) / np.percentile(_lado, 5):.2f})"
)
if sum(v > 2 for v in cuenta.values()) or libres != 2 * N_directo + N_s:
    print("AVISO: malla no conforme, usar más divisiones en el directo (N_directo >= 32)")

# =============================================================================
# 7. GRÁFICO
# =============================================================================
if __name__ == "__main__":
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    colores = np.array(["lightsteelblue", "khaki", "palegreen"])
    ax.add_collection3d(
        Poly3DCollection(nodos[elementos], facecolors=colores[grupo], edgecolors="k", linewidths=0.5, alpha=0.95)
    )
    ax.scatter(*nodos_curva.T, color="g", s=12)
    mn, mx = nodos.min(axis=0), nodos.max(axis=0)
    ax.set_xlim(mn[0], mx[0]); ax.set_ylim(mn[1], mx[1]); ax.set_zlim(mn[2], mx[2])
    ax.set_box_aspect(mx - mn)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    plt.show()
