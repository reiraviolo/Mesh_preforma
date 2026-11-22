"""
ensamblador.py — genera la malla real de cada tramo (tubo recto) y de cada
accesorio (Tee/Codo/Reduccion), ya rotados y trasladados a su posición en
el espacio global según lo que calculó red.calcular_posiciones(), y las
junta en una sola lista de piezas lista para graficar.

OJO — alcance de esta versión: cada pieza se genera y posiciona por
separado; no se fusionan los nodos que coinciden en cada unión en una
malla única (eso es un paso más, de "unir mallas", parecido a unir_mallas()
de tee.py pero para N piezas en vez de 2). Acá el objetivo es que la red
completa se vea bien armada en 3D.

Ejes locales de cada pieza: las cuatro entran por +X (verificado corriendo
cada generador y mirando el resultado), así la sección transversal de
entrada queda siempre en el mismo plano Y-Z para las cuatro, y acoplan
bien cuando más adelante se fusionen los nodos en cada unión:

    tubo.py       : eje de entrada local = +X. Sin desvío (sección circular).
    reduccion.py  : eje de entrada local = +X. Sin desvío.
    codo.py       : eje de entrada local = +X. Se desvía en el plano X-Z,
                    es decir alrededor del eje local Y.
    tee.py        : eje "directo" (tronco) local = +X. La derivación sale
                    en el plano X-Z (mismo eje normal que el codo: Y local).
                    El centro (x=0) es el punto donde tronco y derivación
                    se cruzan — ahí va la posición del nodo.
"""

import numpy as np

from transformaciones import alinear_pieza
from red import calcular_posiciones, _eje_desvio
from tubo import generar_malla_tubo
from codo import generar_malla_codo
from reduccion import generar_malla_reduccion
from tee import generar_malla_tee

# eje local de entrada: +X para las cuatro piezas (ver docstring arriba).
# El eje normal (bisagra del desvío) es -Y, no +Y: tanto codo.py como
# tee.py curvan su tangente de (1,0,0) hacia (cos t, 0, sin t), y esa
# trayectoria corresponde a una rotación de Rodrigues alrededor de -Y (se
# verificó numéricamente comparando la tangente real de salida de la
# malla contra la dirección calculada por red.py: con +Y salía invertida).
EJE_ENTRADA_LOCAL = np.array([1.0, 0.0, 0.0])
EJE_NORMAL_LOCAL = np.array([0.0, -1.0, 0.0])  # igual para codo y tee


def _tangente_codo(params):
    """Distancia (tangente) desde el vértice -la intersección de los ejes
    rectos de entrada y salida, prolongados- hasta el punto donde arranca
    la curva del codo: T = R * tan(angulo/2) (geometría estándar de un
    codo: las dos tangentes desde el vértice al círculo de curvatura
    tienen la misma longitud)."""
    return params["R"] * np.tan(np.deg2rad(params["Angulo"]) / 2)


def _pieza_tramo(tramo, pos_origen, dir_tramo, ref_tramo, N, recorte_inicio=0.0, recorte_fin=0.0):
    """Genera el tubo recto de un tramo, ya rotado y trasladado.

    recorte_inicio/recorte_fin: cuánto se acorta el tramo en cada punta
    porque ahí ya empieza el "cuerpo" propio del accesorio contiguo (su
    malla ya cubre esa distancia, así que el tubo no debe superponerse).
    recorte_inicio además desplaza el punto de partida del tubo esa misma
    distancia hacia adelante (el tubo no arranca en pos_origen sino un
    poco más allá, justo donde termina el accesorio anterior).

    ref_tramo: la "fase" (plano de referencia) heredada de la pieza
    anterior (posiciones[destino]["ref_in"]). Aunque un tubo recto es
    rotacionalmente simétrico en lo geométrico (cualquier giro de su malla
    continua representa el mismo tubo), su malla DISCRETA no lo es: los N
    nodos del perímetro están en ángulos fijos, y para que coincidan
    exactamente (nodo a nodo, no solo "visualmente cerca") con los de la
    pieza contigua -imprescindible para que la malla sirva para FEM- hace
    falta usar la MISMA referencia de fase que esa pieza, no un eje fijo
    cualquiera. Antes se usaba EJE_NORMAL_LOCAL fijo acá, que solo por
    casualidad coincide con la fase correcta en el primerísimo tramo (el
    que sale de la raíz) y en cadenas sin ningún Codo/Tee en el medio; en
    cuanto hay una pieza con Giro != 0 antes, la fase queda desalineada y
    el tramo que sigue no suelda del todo con sus vecinos (aristas libres
    de más, mismo lugar pero sin fusionar)."""
    L = tramo["params"]["L"] - recorte_inicio - recorte_fin
    if L <= 0:
        raise ValueError(
            f"Tramo {tramo['origen']}->{tramo['destino']}: largo {tramo['params']['L']} "
            f"es menor que la suma de los recortes ({recorte_inicio} + {recorte_fin}) "
            f"de los accesorios en sus extremos; revisa L o H_1/H_2/H_3/R/Angulo/L_reduccion."
        )
    OD = tramo["params"]["OD"]
    pos_efectiva = pos_origen + dir_tramo * recorte_inicio
    r = generar_malla_tubo(OD=OD, L=L, N=N, verbose=False)
    nodos_rot = alinear_pieza(
        r["nodos"], EJE_ENTRADA_LOCAL, EJE_NORMAL_LOCAL, dir_tramo, ref_tramo
    )
    nodos_final = nodos_rot + pos_efectiva
    return dict(nodos=nodos_final, elementos=r["elementos"], nombre=f"tramo {tramo['origen']}->{tramo['destino']}")


def _pieza_codo(nodo_id, params, pos, dir_in, ref_in, N):
    OD = params["OD"]  # el OD de un Codo viene del tramo entrante (ver ensamblar())
    r = generar_malla_codo(OD=OD, R=params["R"], angulo=params["Angulo"], N=N, verbose=False)
    # "Giro" se mide respecto del plano XY de la pieza anterior (ref_in),
    # no de una horizontal global (ver docstring de red.py).
    eje_desvio = _eje_desvio(ref_in, dir_in, params.get("Giro", 0))
    nodos_rot = alinear_pieza(r["nodos"], EJE_ENTRADA_LOCAL, EJE_NORMAL_LOCAL, dir_in, eje_desvio)
    nodos_final = nodos_rot + pos
    return dict(nodos=nodos_final, elementos=r["elementos"], nombre=f"codo {nodo_id}"), None


def _pieza_reduccion(nodo_id, params, pos, dir_in, ref_in, N):
    r = generar_malla_reduccion(OD_1=params["OD1"], OD_2=params["OD2"], L=params["L"], N=N, verbose=False)
    # misma fase heredada que el tramo recto (ver nota en _pieza_tramo):
    # la reducción no cambia de dirección, así que su propio ref_in debe
    # pasar sin tocar a lo que venga después, pero también define la fase
    # con la que hay que alinearla a lo que viene ANTES.
    nodos_rot = alinear_pieza(
        r["nodos"], EJE_ENTRADA_LOCAL, EJE_NORMAL_LOCAL, dir_in, ref_in
    )
    nodos_final = nodos_rot + pos
    return dict(nodos=nodos_final, elementos=r["elementos"], nombre=f"reduccion {nodo_id}"), None


def _pieza_tee(nodo_id, params, pos, dir_in, ref_in, N):
    OD = params.get("OD", None)
    # OD_1 (tronco) viene del tramo entrante, igual que en el Codo.
    r = generar_malla_tee(
        OD_1=OD,
        OD_2=params["OD2"],
        H_1=params["H1"],
        H_2=params["H2"],
        H_3=params["H3"],
        Alpha_2=params["Alpha"],
        N_directo=N,
        verbose=False,
    )
    # "Giro" se mide respecto del plano XY de la pieza anterior (ref_in),
    # no de una horizontal global (ver docstring de red.py).
    # tee.py puede, en silencio (solo un print), AMPLIAR H_1/H_3 si la corona
    # de la derivación no entra en el H declarado por el usuario (ver AVISO
    # en parts/tee.py). Si eso pasa, el tramo recto contiguo -que ensamblador
    # recorta usando el H1/H3 NOMINAL, el que escribió el usuario- queda
    # terminando antes (o despues) de donde realmente termina la malla de
    # la Tee: un desacople de POSICION (no de fase) entre ambas piezas que
    # rompe la coincidencia exacta de nodos exigida para FEM. En vez de
    # tolerarlo en silencio, se lo transforma en un error explicito y
    # accionable: que el usuario suba H_1/H_3 al minimo real.
    H1_real, H3_real = r["divisiones"]["H_1"], r["divisiones"]["H_3"]
    if H1_real != params["H1"] or H3_real != params["H3"]:
        raise ValueError(
            f"Nodo {nodo_id} (Tee): la corona de la derivacion no entra en "
            f"H1={params['H1']:g}/H3={params['H3']:g}; para que la malla quede "
            f"exactamente conforme con los tramos contiguos (requisito de FEM) "
            f"hace falta declarar H1>={H1_real:g} y H3>={H3_real:g}."
        )
    eje_desvio = _eje_desvio(ref_in, dir_in, params.get("Giro", 0))
    nodos_rot = alinear_pieza(r["nodos"], EJE_ENTRADA_LOCAL, EJE_NORMAL_LOCAL, dir_in, eje_desvio)
    nodos_final = nodos_rot + pos
    N_rama = r["divisiones"]["N_s"]  # el que tee.py elige para la derivación
    return dict(nodos=nodos_final, elementos=r["elementos"], nombre=f"tee {nodo_id}"), N_rama


def ensamblar(nodos, tramos, N, direccion_inicial=(1.0, 0.0, 0.0), verbose=True):
    """
    Genera todas las piezas (tramos + accesorios) de la red, posicionadas
    y orientadas en el espacio global.

    N : número de divisiones en el perímetro a usar en el tronco principal
        y en cada tubo/codo/reducción. En la derivación de una Tee, ese
        número se recalcula solo (el N_s que tee.py elige para que los
        elementos queden lo más cuadrados posible) y es ESE número el que
        se usa de ahí en adelante en toda la rama de la derivación —
        tramo siguiente y cualquier accesorio más adelante en esa rama.

    Devuelve una lista de piezas, cada una {"nodos", "elementos", "nombre"}.
    """
    posiciones = calcular_posiciones(nodos, tramos, direccion_inicial)

    salientes = {nid: [] for nid in nodos}
    con_entrante = set()
    for t in tramos:
        salientes[t["origen"]].append(t)
        con_entrante.add(t["destino"])
    raices = [nid for nid in nodos if nid not in con_entrante]
    raiz = raices[0]

    piezas = []
    N_in = {raiz: N}  # N con el que se genera lo que ENTRA a cada nodo
    OD_in = {}  # OD del tramo que entra a cada nodo (Codo/Tee lo necesitan
    #              y no viene en los parametros del propio nodo)
    pendientes = [raiz]

    while pendientes:
        actual = pendientes.pop()
        n_actual = N_in[actual]
        tipo = nodos[actual]["tipo"]
        params = dict(nodos[actual]["params"])  # copia, le agregamos OD si hace falta
        if actual in OD_in:
            params["OD"] = OD_in[actual]
        pos = posiciones[actual]["pos"]
        dir_in = posiciones[actual]["dir_in"]
        ref_in = posiciones[actual]["ref_in"]

        N_rama = n_actual  # por defecto Codo/Reduccion no cambian el N
        # Fase real (ya rotada "Giro") de TODA la cara de salida de "actual",
        # si es un Codo/Tee: es una sola pieza rigida, asi que esta misma
        # fase es la que ve cualquier tramo que salga de aca (tronco o
        # derivacion da igual, ver nota mas abajo). None si "actual" no
        # impone ninguna fase propia (Libre, Reduccion: dejan pasar la que
        # traigan, sin fijar nada).
        eje_desvio_actual = None

        if tipo == "Codo" and dir_in is not None:
            # el codo se ubica en el vertice de la interseccion de los ejes
            # rectos de entrada y salida; su propia malla, en cambio,
            # arranca en el punto de tangencia, retrocedido "T" desde el
            # vertice en -dir_in.
            T = _tangente_codo(params)
            pieza, _ = _pieza_codo(actual, params, pos - dir_in * T, dir_in, ref_in, n_actual)
            piezas.append(pieza)
            eje_desvio_actual = _eje_desvio(ref_in, dir_in, params.get("Giro", 0))
        elif tipo == "Reduccion" and dir_in is not None:
            # el punto de inserción de la reduccion es su propio punto de
            # entrada (pos ya es ese punto, sin corrimiento).
            pieza, _ = _pieza_reduccion(actual, params, pos, dir_in, ref_in, n_actual)
            piezas.append(pieza)
        elif tipo == "Tee" and dir_in is not None:
            # igual que la reduccion: pos ya es el punto donde se cruzan
            # el tronco y la derivacion, que es el origen local (x=0) de
            # tee.py, sin corrimiento.
            pieza, N_rama_tee = _pieza_tee(actual, params, pos, dir_in, ref_in, n_actual)
            piezas.append(pieza)
            N_rama = N_rama_tee
            eje_desvio_actual = _eje_desvio(ref_in, dir_in, params.get("Giro", 0))
            if verbose:
                print(f"[ensamblador] Tee {actual}: tronco N={n_actual} -> derivacion N_s={N_rama}")
        # tipo == "Libre": no genera pieza propia, solo deja pasar N_in

        for t in salientes[actual]:
            dir_tramo = posiciones[t["destino"]]["dir_in"]
            # Fase con la que hay que alinear ESTE tramo para que coincida
            # nodo a nodo con la cara de entrada del accesorio al que
            # llega (ver nota en _pieza_tramo). Para Codo/Tee esa cara de
            # entrada no usa el "ref_in" tal cual: usa "eje_desvio" (ref_in
            # ya rotado "Giro" grados), que es lo mismo que usa
            # _pieza_codo()/_pieza_tee() como normal objetivo de TODA la
            # pieza (entrada y desvío son una sola rotación rígida, no se
            # pueden separar). Para Reduccion/Libre no hay Giro que rote
            # nada, así que el "ref_in" sin modificar es la fase correcta
            # (coincide con lo que usa _pieza_reduccion()).
            tipo_destino_fase = nodos[t["destino"]]["tipo"]
            ref_in_destino = posiciones[t["destino"]]["ref_in"]
            dir_in_destino = posiciones[t["destino"]]["dir_in"]
            if tipo_destino_fase in ("Codo", "Tee"):
                giro_destino = nodos[t["destino"]]["params"].get("Giro", 0)
                ref_tramo = _eje_desvio(ref_in_destino, dir_in_destino, giro_destino)
            else:
                ref_tramo = ref_in_destino
            es_derivacion = tipo == "Tee" and t["params"].get("D")
            n_tramo = N_rama if es_derivacion else n_actual

            # Si TANTO el origen como el destino de este tramo son, cada
            # uno, un Codo/Tee (un accesorio enchufado directo a otro, sin
            # ninguna Reduccion/Libre de por medio), ambos extremos imponen
            # su propia fase (ver eje_desvio_actual / ref_tramo arriba) y
            # el tramo -recto, con una sola fase uniforme en toda su
            # longitud- solo puede calzar exacto con los dos a la vez si
            # esas dos fases forzadas difieren en un múltiplo entero del
            # paso angular de la malla (360°/n_tramo); si no, no existe
            # ninguna fase que sirva para ambos lados a la vez y los nodos
            # quedan desalineados en una de las dos uniones (requisito de
            # FEM incumplido). Esto no es un bug de signo o de que falte
            # propagar algo: es geometria -un caño recto no puede "destorcer"
            # un desfasaje que no es multiplo del paso de grilla- así que se
            # avisa con un error accionable en vez de fusionar mal en
            # silencio.
            if eje_desvio_actual is not None and tipo_destino_fase in ("Codo", "Tee"):
                paso = 360.0 / n_tramo
                cos_ang = np.clip(np.dot(eje_desvio_actual, ref_tramo)
                                   / (np.linalg.norm(eje_desvio_actual) * np.linalg.norm(ref_tramo)), -1.0, 1.0)
                desfase = np.degrees(np.arccos(cos_ang))
                # el sentido (signo) no importa para esto, solo que el
                # valor absoluto sea multiplo de "paso".
                resto = desfase % paso
                resto = min(resto, paso - resto)
                if resto > 0.05:
                    raise ValueError(
                        f"Nodo {actual} ({tipo}) -> Nodo {t['destino']} ({tipo_destino_fase}): "
                        f"un accesorio enchufado directo a otro (sin Reduccion/Libre de por "
                        f"medio) necesita que la malla calce exacta en AMBAS caras del tramo "
                        f"recto que los une, pero las fases reales de esas dos caras difieren "
                        f"{desfase:.2f}°, que no es multiplo de los {paso:.2f}° de paso de la "
                        f"malla (N={n_tramo}). Para que los nodos coincidan exactos (FEM), "
                        f"ajusta el 'Giro' de uno de los dos accesorios, o N, de forma que la "
                        f"diferencia entre sus Giros sea multiplo de {paso:.2f}°, o intercala un "
                        f"tramo con una Reduccion/nodo Libre entre ambos."
                    )

            # recorte al inicio del tramo: cuanto ya cubre, en esta punta,
            # el accesorio del que sale ("actual")
            if tipo == "Tee":
                recorte_inicio = params["H2"] if es_derivacion else params["H3"]
            elif tipo == "Codo":
                recorte_inicio = _tangente_codo(params)
            elif tipo == "Reduccion":
                recorte_inicio = params["L"]
            else:  # Libre
                recorte_inicio = 0.0

            # recorte al final del tramo: cuanto va a cubrir, en esta otra
            # punta, el accesorio al que entra (si tiene uno)
            tipo_destino = nodos[t["destino"]]["tipo"]
            params_destino = nodos[t["destino"]]["params"]
            if tipo_destino == "Tee":
                recorte_fin = params_destino["H1"]
            elif tipo_destino == "Codo":
                recorte_fin = _tangente_codo(params_destino)
            else:  # Reduccion (sin recorte en su entrada) o Libre
                recorte_fin = 0.0

            piezas.append(
                _pieza_tramo(t, pos, dir_tramo, ref_tramo, n_tramo, recorte_inicio, recorte_fin)
            )
            N_in[t["destino"]] = n_tramo
            OD_in[t["destino"]] = t["params"]["OD"]
            pendientes.append(t["destino"])

    return piezas


def _agrupar_soldadura(elementos, es_soldadura):
    """Agrupa los nodos de soldadura en zonas (una zona = un anillo de
    soldadura entre dos piezas contiguas), uniendo los que quedan conectados
    por una arista de algún elemento de la malla fusionada.

    elementos     : (E,4) int, conectividad de la malla ya fusionada
    es_soldadura  : (N,) bool, True en los nodos que vinieron de fusionar
                    >=2 nodos originales (es decir, nodos de unión)

    Devuelve una lista de arrays de ids de nodo, uno por zona, ordenadas
    por su primer id (para que la salida sea reproducible).
    """
    padre = {i: i for i in range(len(es_soldadura)) if es_soldadura[i]}

    def raiz(i):
        while padre[i] != i:
            padre[i] = padre[padre[i]]
            i = padre[i]
        return i

    def unir(i, j):
        ri, rj = raiz(i), raiz(j)
        if ri != rj:
            padre[ri] = rj

    # aristas de cada cuadrilatero (o triangulo [a,b,c,c]); si ambos
    # extremos son nodos de soldadura, pertenecen a la misma zona
    for a, b, c, d in elementos:
        for u, v in ((a, b), (b, c), (c, d), (d, a)):
            if u != v and es_soldadura[u] and es_soldadura[v]:
                unir(u, v)

    zonas = {}
    for i in padre:
        zonas.setdefault(raiz(i), []).append(i)

    zonas_ordenadas = sorted(zonas.values(), key=min)
    return [np.array(sorted(z), dtype=int) for z in zonas_ordenadas]


def fusionar_piezas(piezas, tol=1e-3):
    """Junta todas las piezas en una sola malla, fusionando los nodos que
    coinciden (dentro de la tolerancia `tol`) en cada unión entre piezas
    contiguas — así la malla final queda con nodos compartidos ahí, en vez
    de piezas sueltas que solo coinciden visualmente.

    Esos nodos fusionados son, físicamente, los nodos de soldadura: el
    punto donde dos piezas (tubo/codo/reducción/tee) se sueldan entre sí.

    tol : distancia máxima entre dos nodos para considerarlos "el mismo"
          (mismas unidades que las coordenadas de la malla, p. ej. mm).

    Devuelve (nodos, elementos, zonas_soldadura):
      nodos             : (N,3) de la malla fusionada completa
      elementos         : (E,4) de la malla fusionada completa
      zonas_soldadura   : lista de arrays de ids de nodo; cada array es
                          un anillo de soldadura (una unión entre piezas)
    """
    nodos_todos = []
    elementos_todos = []
    offset = 0
    for p in piezas:
        n = len(p["nodos"])
        nodos_todos.append(p["nodos"])
        elementos_todos.append(np.asarray(p["elementos"], dtype=int) + offset)
        offset += n

    nodos = np.vstack(nodos_todos)
    elementos = np.vstack(elementos_todos)

    # agrupa nodos "iguales" redondeando sus coordenadas a la tolerancia
    # dada; np.unique con axis=0 asigna el mismo id a cada grupo repetido.
    clave = np.round(nodos / tol).astype(np.int64)
    _, idx_unico, inverso = np.unique(clave, axis=0, return_index=True, return_inverse=True)

    nodos_fusionados = nodos[idx_unico]
    elementos_fusionados = inverso[elementos]

    # un nodo fusionado es "de soldadura" si recibió 2 o más nodos
    # originales (es decir, nodos de piezas distintas que coincidieron ahí)
    conteo = np.bincount(inverso, minlength=len(nodos_fusionados))
    es_soldadura = conteo > 1
    zonas_soldadura = _agrupar_soldadura(elementos_fusionados, es_soldadura)

    return nodos_fusionados, elementos_fusionados, zonas_soldadura


if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt

    sys.path.insert(0, ".")
    from mermaid_reader import leer_mermaid
    import graficador

    if len(sys.argv) < 2:
        print("Uso: python ensamblador.py archivo.md [N]")
        sys.exit(1)

    N = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    nodos, tramos = leer_mermaid(sys.argv[1])
    piezas = ensamblar(nodos, tramos, N=N)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    for p in piezas:
        graficador.graficar_malla(p["nodos"], p["elementos"], ax=ax, mostrar=False)
    plt.title(f"Red ensamblada ({len(piezas)} piezas, N={N})")
    plt.show()
