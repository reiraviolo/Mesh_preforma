"""Malla preformada de una tee por BLOQUES (solo cuadriláteros).

Idea (ver croquis): las superficies se dividen en bloques estructurados.

  * Tubo directo: anillos con N_directo nodos equiespaciados en la circunferencia
    (las dos bocas quedan equidistanciadas). La posición de los anillos a lo largo
    del eje SE ACOMODA para que haya una "ventana" rectangular alrededor de la
    derivación: x = xa ... xb y dos generatrices a cada lado.
  * Ventana: se rellena con una corona de K capas (malla en O) entre el borde de la
    ventana y la curva de intersección.
  * Derivación: N_salida nodos equiespaciados en la circunferencia (boca equidistanciada).
    Sus generatrices arrancan EXACTAMENTE en los nodos de la curva, que son los mismos
    nodos de la corona del directo (coincidentes por construcción).

Pasos:  1. Parámetros  2. Curva de intersección  3. Tubo directo (grilla + ventana)
        4. Corona  5. Derivación  6. Control  7. Gráfico
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
H_3 = 200  # de intersección a salida final
Alpha_2 = 90  # ángulo de la derivación [grados]

N_directo = 70  # divisiones de la boca del directo (se redondea a múltiplo de 4)
CAPAS = 1  # capas de elementos entre la ventana y la curva (corona)
SUAVIZADO = 60  # iteraciones de suavizado: la grilla se deforma para adaptarse a la curva
DIST_MIN = 1.0  # los nodos del borde de la ventana no se acercan a la curva a menos de DIST_MIN·tamaño
BANDA = 3  # anillos a cada lado de la ventana que pueden deformarse
TRANSICION = 3.0  # largo de la zona de transición de tamaños junto a la ventana (en tamaños de ventana)
ALARGAMIENTO_MAX = 3.0  # máx. relación entre el tamaño en x de la ventana y el del directo
N_X = None  # celdas de la ventana en x (None = las mínimas: N_directo/10, así hay pocos nodos agregados)
LIMITE_NODOS = 1.15  # la derivación no puede tener más de este factor de divisiones que el directo

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
# 3.0 N_directo mínimo para que la ventana respete LIMITE_NODOS (si hace falta se sube
#     de a 4 y se avisa)
def ventana(N):
    dth_, tam_ = 2 * np.pi / N, R1 * 2 * np.pi / N
    n_th_ = 2 * int(np.ceil((np.pi / 2 - th_min) / dth_ + CAPAS))
    ancho_ = x_max - x_min + 2 * CAPAS * tam_
    n_x_nat = int(np.ceil(ancho_ / tam_))
    n_x_ = min(n_x_nat, int(LIMITE_NODOS * N / 2) - n_th_, max(4, 2 * int(np.ceil(N / 20))))
    n_x_ = max(n_x_, int(np.ceil(ancho_ / (ALARGAMIENTO_MAX * tam_))))  # celdas no demasiado alargadas
    n_x_ += n_x_ % 2
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

# 3.2 anillos. La ventana tiene n_x celdas (par) y cubre la curva + margen. Como la
#     derivación tiene N_s = 2 (n_th + n_x) divisiones, n_x se limita para que
#     N_s <= LIMITE_NODOS * N_directo; si hace falta, las celdas de la ventana se
#     alargan en x. Fuera de la ventana el tamaño crece/decrece gradualmente.
margen = CAPAS * tam
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


X_izq = anillos_con_rampa(-H_1, xa, tam, tam_x, True)
X_der = anillos_con_rampa(xb, H_3, tam_x, tam, False)
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
        if not dentro:
            cuads_grilla.append(
                [id_grilla(i, j), id_grilla(i, j + 1), id_grilla(i + 1, j + 1), id_grilla(i + 1, j)]
            )
cuads_grilla = np.array(cuads_grilla)

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
N_s = len(lazo)  # divisiones de la derivación
lazo = lazo[n_x // 2 :] + lazo[: n_x // 2]  # ahora el nodo 0 corresponde a φ = 0
assert N_s == 2 * (n_x + n_th)

# Ángulos φ de los nodos sobre la curva. Las 4 esquinas de la ventana se unen con los
# puntos de la curva que están en la misma dirección desde el centro de la ventana; entre
# esquinas los nodos se reparten por longitud de arco (en el desarrollo x, R·θ), así los
# elementos de la corona quedan lo menos deformados posible.
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
p_esq = [n_x // 2, n_x // 2 + n_th, N_s // 2, 3 * n_x // 2 + n_th, 3 * n_x // 2 + 2 * n_th, N_s, N_s + n_x // 2]
phi_c = np.zeros(N_s)
for k in range(6):
    a0 = np.interp(phi_esq[k], phi_ext, arco)
    a1 = np.interp(phi_esq[k + 1], phi_ext, arco)
    for p in range(p_esq[k], p_esq[k + 1]):
        phi_c[p % N_s] = np.interp(a0 + (a1 - a0) * (p - p_esq[k]) / (p_esq[k + 1] - p_esq[k]), arco, phi_ext) % (2 * np.pi)
nodos_curva = punto_curva(phi_c)

# ángulos φ de la boca de la derivación: equiespaciados, alineados con los de la curva
phi_boca = np.arange(N_s) * 2 * np.pi / N_s
desfase = np.angle(np.mean(np.exp(1j * (phi_c - phi_boca))))
phi_boca = phi_boca + desfase

id_grilla_borde = np.array([id_grilla(i, j) for i, j in lazo])
n_g = len(nodos_grilla)
n_corona_int = N_s * (CAPAS - 1)  # nodos intermedios de la corona
id_curva = n_g + n_corona_int + np.arange(N_s)  # numeración de los nodos de la curva


def a_plano(p):
    """(x, θ) de puntos sobre el cilindro directo."""
    return np.stack([p[..., 0], np.arctan2(p[..., 2], p[..., 1])], axis=-1)


a_borde = np.array([[X[i], j * dth] for i, j in lazo])  # (x, θ) con θ continuo (j puede salirse de 0..N-1)
a_curva = a_plano(nodos_curva)
nodos_int = []
capa_id = np.zeros((N_s, CAPAS + 1), dtype=int)
capa_id[:, 0] = id_grilla_borde
capa_id[:, CAPAS] = id_curva
for k in range(1, CAPAS):
    s = k / CAPAS
    plano = (1 - s) * a_borde + s * a_curva  # interpolación en el desarrollo (x, θ)
    nodos_int.append(
        np.stack([plano[:, 0], R1 * np.cos(plano[:, 1]), R1 * np.sin(plano[:, 1])], axis=1)
    )
    capa_id[:, k] = n_g + (k - 1) * N_s + np.arange(N_s)
nodos_int = np.vstack(nodos_int) if nodos_int else np.zeros((0, 3))
cuads_corona = np.array(
    [
        [capa_id[p, k], capa_id[(p + 1) % N_s, k], capa_id[(p + 1) % N_s, k + 1], capa_id[p, k + 1]]
        for p in range(N_s)
        for k in range(CAPAS)
    ]
)

# =============================================================================
# 5. DERIVACIÓN: generatrices desde los nodos de la curva hasta la boca (t = H_2)
# =============================================================================
t_curva = nodos_curva @ eje_s  # posición axial de cada nodo de la curva
tam_s = 2 * np.pi * R2 / N_s  # tamaño circunferencial de la derivación
tam_axial = tam_s  # elementos cuadrados en la derivación (mismo tamaño axial que circunferencial)
M = max(1, int(round(np.mean(H_2 - t_curva) / tam_axial)))  # capas a lo largo del eje
nodos_der = []
der_id = np.zeros((N_s, M + 1), dtype=int)
der_id[:, 0] = id_curva
base = n_g + n_corona_int + N_s
dphi = np.angle(np.exp(1j * (phi_boca - phi_c)))  # giro de las generatrices hasta la boca
for m in range(1, M + 1):
    u = m / M
    phi = phi_c + dphi * u  # la última capa (boca) queda equiespaciada y plana en t = H_2
    t = t_curva + (H_2 - t_curva) * u
    nodos_der.append(
        t[:, None] * eje_s + (R2 * np.cos(phi))[:, None] * np.array([0, 1.0, 0]) + (R2 * np.sin(phi))[:, None] * perp_s
    )
    der_id[:, m] = base + (m - 1) * N_s + np.arange(N_s)
nodos_der = np.vstack(nodos_der)
cuads_der = np.array(
    [
        [der_id[j, m], der_id[(j + 1) % N_s, m], der_id[(j + 1) % N_s, m + 1], der_id[j, m + 1]]
        for j in range(N_s)
        for m in range(M)
    ]
)

# ---- ensamble --------------------------------------------------------------
nodos = np.vstack([nodos_grilla, nodos_int, nodos_curva, nodos_der])
elementos = np.vstack([cuads_grilla, cuads_corona, cuads_der])
grupo = np.concatenate(
    [np.zeros(len(cuads_grilla), int), np.ones(len(cuads_corona), int), np.full(len(cuads_der), 2)]
)


# ---- suavizado: la grilla del directo se deforma para adaptarse a la curva -------------
def suavizar(nodos, cuads, libres_ids, iteraciones, w=0.6, cuidar=(), curva=None, dist_min=0.0):
    """Laplaciano sobre el cilindro directo; los nodos de `libres_ids` se mueven, el resto no."""
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
            q = nodos[i] + w * (m - nodos[i])
            r = np.hypot(q[1], q[2])
            q = np.array([q[0], q[1] * R1 / r, q[2] * R1 / r])  # de nuevo sobre el cilindro
            if i in cuidar and np.linalg.norm(curva - q, axis=1).min() < dist_min:
                continue  # un nodo del borde de la ventana no se acerca demasiado a la curva
            nuevo[i] = q
        nodos = nuevo
    return nodos


cuads_directo = np.vstack([cuads_grilla, cuads_corona])
anillo_de = np.arange(n_g) // N_directo
libres_ids = [
    i for i in range(n_g) if i1 - BANDA <= anillo_de[i] <= i2 + BANDA and 0 < anillo_de[i] < n_anillos - 1
]
borde_ventana = set(id_grilla_borde.tolist())
libres_ids = libres_ids + list(range(n_g, n_g + n_corona_int))  # capas intermedias de la corona
if SUAVIZADO:
    nodos = suavizar(
        nodos, cuads_directo, libres_ids, SUAVIZADO, cuidar=borde_ventana, curva=nodos_curva, dist_min=DIST_MIN * tam
    )

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

if N_s > LIMITE_NODOS * N_directo:
    print(
        f"AVISO: la derivación tiene {N_s / N_directo - 1:.0%} más divisiones que el directo (límite "
        f"{LIMITE_NODOS - 1:.0%}); con esta geometría no se puede bajar sin deformar demasiado los elementos"
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
