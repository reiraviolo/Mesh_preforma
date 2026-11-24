"""
red.py — resuelve la geometría global de la red a partir de lo que devuelve
mermaid_reader.leer_mermaid(): para cada nodo calcula su posición 3D y la
dirección del tramo que entra, recorriendo el grafo desde la raíz.

No genera ninguna malla todavía: es el paso intermedio para poder verificar
(con prints o graficando los ejes) que el ruteo da lo esperado antes de
invertir tiempo en mallar.

Convención de direcciones:
- El primer tramo (el que sale del nodo raíz, sin tramo entrante) usa una
  dirección inicial fija por convención: (1, 0, 0).
- Codo: gira la dirección entrante por "Angulo" grados, alrededor de un eje
  perpendicular calculado con "Giro" grados respecto del plano XY de la
  PIEZA ANTERIOR (no de una horizontal global fija: ver "ref" más abajo).
- Tee: el tramo marcado con el flag S sigue derecho (misma dirección que
  entra); el marcado con el flag D sale desviado "Alpha" grados, con el
  mismo esquema de "Giro" que el Codo.
- Reduccion: sigue derecho (solo cambia el diámetro, no la dirección).
- Libre: nodo hoja (no tiene salientes) o raíz (no tiene entrante).

Convención de "Giro" (importante):
- "Giro" no se mide contra una horizontal absoluta del espacio: se mide
  respecto del plano XY de la pieza anterior (su propio plano de desvío,
  o -si la pieza anterior era recta- el que traía de más atrás). Por eso
  cada nodo lleva, ademas de "dir_in" (direccion del tramo entrante), un
  "ref_in": el vector perpendicular a dir_in que define donde esta el
  Giro=0 para el accesorio de ESE nodo, heredado de la pieza anterior.
- Una pieza recta (Reduccion, tramo derecho, salida "S" de una Tee) no
  cambia la direccion, asi que el ref_in pasa sin modificar.
- Un Codo o la derivacion "D" de una Tee SI desvian la direccion: el nuevo
  ref que se propaga aguas abajo es el propio eje de desvio de esa pieza
  (la normal a su plano de giro), que por construccion queda perpendicular
  tanto a la direccion entrante como a la saliente.
- Solo el primerisimo accesorio de la red (el que sigue directo a la raiz,
  sin ninguna pieza antes) no tiene de donde heredar un ref: para ese caso
  se arma uno por convencion con _ref_inicial() (perpendicular a la
  direccion inicial, usando el eje Z global como referencia).
"""

import numpy as np

from transformaciones import rotar_eje

GRADOS_ESPERADOS = {"Tee": 3, "Codo": 2, "Reduccion": 2, "Libre": 1}


def _validar_grados(nodos, tramos):
    grado = {nid: 0 for nid in nodos}
    for t in tramos:
        grado[t["origen"]] += 1
        grado[t["destino"]] += 1

    errores = []
    for nid, info in nodos.items():
        esperado = GRADOS_ESPERADOS.get(info["tipo"])
        if esperado is not None and grado.get(nid, 0) != esperado:
            errores.append(
                f"  - Nodo {nid} ({info['tipo']}): tiene {grado.get(nid, 0)} "
                f"tramo(s), se esperaba(n) {esperado}."
            )
    if errores:
        raise ValueError("Topologia inconsistente:\n" + "\n".join(errores))


def _ref_inicial(direccion):
    """Arma, por convención, el "ref" (bisagra de referencia para Giro=0)
    de la primerísima pieza de la red: la que sale directo de la raíz, sin
    ninguna pieza antes de la que heredar una orientación.

    "ref" es el EJE (bisagra) alrededor del cual se gira `direccion` para
    obtener la salida del Codo/Tee, no la dirección de salida en sí. Por
    convención, Giro=0 usa como bisagra el eje Z global (proyectado
    perpendicular a `direccion`, por si `direccion` no es exactamente
    horizontal): una bisagra vertical da un codo que dobla en el plano
    horizontal (sale hacia +Y si `direccion` es +X), que es el default
    esperado. Si `direccion` ya es ~vertical, Z no sirve de referencia (casi
    paralelo), así que ahí se usa X global en su lugar.
    """
    direccion = direccion / np.linalg.norm(direccion)
    z = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(direccion, z)) > 0.999:
        ref = np.array([1.0, 0.0, 0.0])  # direccion casi paralela a Z: usar X
        ref = ref - np.dot(ref, direccion) * direccion
    else:
        ref = z - np.dot(z, direccion) * direccion  # Z proyectado perpendicular a direccion
    return ref / np.linalg.norm(ref)


def _eje_desvio(ref_in, dir_in, giro_grados):
    """Calcula el eje alrededor del cual hay que girar dir_in para obtener
    la salida de un Codo o la derivación de una Tee.

    `ref_in` es el plano XY de la pieza anterior (ver docstring del módulo):
    el vector perpendicular a dir_in que define dónde está Giro=0 para ESTE
    accesorio, heredado del nodo anterior (o armado por _ref_inicial() si es
    el primer accesorio de toda la red).

    Simplemente rota `ref_in` `giro_grados` grados ALREDEDOR de dir_in (no
    alrededor de Y ni de Z) usando rotar_eje(). El resultado es el eje
    perpendicular real sobre el que después se gira dir_in por Angulo o
    Alpha, en _direccion_salida().

    Convención de sentido: con `ref_in` = bisagra vertical (Z) por defecto
    (ver _ref_inicial), Giro=0 da una salida horizontal y Giro=90 una
    salida vertical (+Z) — verificado numéricamente. Esto queda con la
    rotación de rotar_eje() SIN negar el ángulo (regla de la mano derecha
    tal cual, alrededor de dir_in).
    """
    dir_in = dir_in / np.linalg.norm(dir_in)
    return rotar_eje(ref_in, dir_in, giro_grados)


def _direccion_salida(tipo_nodo, dir_in, ref_in, params_nodo, params_tramo):
    """Calcula hacia donde apunta un tramo que SALE de un nodo, a partir de
    la direccion del tramo que ENTRA a ese mismo nodo (dir_in) y de su
    "ref_in" (plano XY heredado de la pieza anterior, ver docstring del
    módulo).

    Devuelve (dir_salida, ref_salida): ref_salida es lo que se propaga como
    ref_in al siguiente nodo (sin cambios si la pieza es recta; el propio
    eje de desvío si la pieza desvía la dirección).

    Cada tipo de accesorio tiene su propia regla:
    """
    if tipo_nodo == "Reduccion":
        # Solo cambia el diametro, la tuberia sigue derecho: misma direccion
        # y mismo plano de referencia (no hay desvio que lo rote).
        return dir_in, ref_in

    if tipo_nodo == "Codo":
        # Un solo tramo de salida: se desvia "Angulo" grados respecto de la
        # entrada, en el plano que define "Giro" (medido desde ref_in, el
        # plano XY de la pieza anterior). El eje de desvio resultante queda
        # perpendicular tanto a dir_in como a dir_salida, asi que es el
        # ref natural para seguir propagando aguas abajo.
        angulo = params_nodo["Angulo"]
        giro = params_nodo.get("Giro", 0)
        eje = _eje_desvio(ref_in, dir_in, giro)
        return rotar_eje(dir_in, eje, angulo), eje

    if tipo_nodo == "Tee":
        # Dos tramos de salida posibles: el marcado "S" sigue derecho
        # (misma direccion que entra), el marcado "D" es la derivacion, que
        # se desvia "Alpha" grados en el plano que define "Giro" (mismo
        # esquema que el Codo, pero con Alpha en vez de Angulo, y midiendo
        # Giro desde el mismo ref_in del tronco).
        #
        # El "ref" que se propaga en AMBOS casos es el mismo "eje" (ref_in
        # ya rotado "Giro" grados): aunque el tronco "S" no cambia de
        # direccion, la pieza de la Tee es UN SOLO cuerpo rigido (tronco +
        # derivacion generados y rotados juntos en ensamblador._pieza_tee,
        # con esa misma rotacion determinada por "eje"), asi que la cara de
        # salida del tronco queda, en los hechos, orientada segun "eje" y
        # no segun el "ref_in" sin rotar -como el propio "Giro" de la Tee
        # ya mueve la fase real de esa cara. Propagar el ref_in sin rotar
        # haria que el tramo recto que sigue (generado con ese ref como
        # fase) no coincida nodo a nodo con la cara real del tronco.
        giro = params_nodo.get("Giro", 0)
        eje = _eje_desvio(ref_in, dir_in, giro)
        if params_tramo.get("S"):
            return dir_in, eje
        if params_tramo.get("D"):
            alpha = params_nodo["Alpha"]
            return rotar_eje(dir_in, eje, alpha), eje
        raise ValueError(
            "Tramo saliente de una Tee sin flag 'S' (directo) ni 'D' "
            "(derivacion); no se puede saber que direccion darle."
        )

    raise ValueError(f"No se sabe como calcular salidas para el tipo '{tipo_nodo}'")


def calcular_posiciones(nodos, tramos, direccion_inicial=(1.0, 0.0, 0.0)):
    """
    Recorre el grafo (debe ser un arbol: sin loops, una sola raiz) y calcula
    la posicion 3D de cada nodo junto con la direccion del tramo entrante.

    Parametros
    ----------
    nodos, tramos : lo que devuelve mermaid_reader.leer_mermaid().
    direccion_inicial : direccion (unitaria o no) del primer tramo, por
        convencion ya que no hay ningun tramo anterior del que partir.

    Devuelve
    --------
    dict {id: {"pos": np.array(3,), "dir_in": np.array(3,) o None,
                "ref_in": np.array(3,) o None}}
    (dir_in/ref_in son None solo para la raiz). "ref_in" es el plano XY de
    la pieza anterior: lo que usa el Codo/Tee de ESE nodo para medir su
    "Giro" (ver docstring del módulo).
    """
    # 1) Chequeo de topología: cada nodo tiene la cantidad de tramos que su
    #    tipo de accesorio permite (Tee=3, Codo=2, Reduccion=2, Libre=1).
    _validar_grados(nodos, tramos)

    # 2) Listas de adyacencia: para cada nodo, sus tramos salientes; y el
    #    conjunto de nodos que SÍ tienen un tramo entrante (todos menos la
    #    raíz, que es justamente el único nodo sin tramo entrante).
    salientes = {nid: [] for nid in nodos}
    con_entrante = set()
    for t in tramos:
        salientes[t["origen"]].append(t)
        con_entrante.add(t["destino"])

    raices = [nid for nid in nodos if nid not in con_entrante]
    if len(raices) != 1:
        raise ValueError(
            f"Se esperaba exactamente una raiz (nodo sin tramo entrante), "
            f"se encontraron {len(raices)}: {raices}"
        )
    raiz = raices[0]

    # 3) Recorrido del árbol (DFS con pila) empezando por la raíz, en
    #    posición (0,0,0) y sin dirección entrante todavía (None).
    resultado = {raiz: {"pos": np.zeros(3), "dir_in": None, "ref_in": None}}
    pendientes = [raiz]

    while pendientes:
        actual = pendientes.pop()
        pos_actual = resultado[actual]["pos"]
        dir_in_actual = resultado[actual]["dir_in"]
        ref_in_actual = resultado[actual]["ref_in"]
        tipo_actual = nodos[actual]["tipo"]
        params_actual = nodos[actual]["params"]

        # Por cada tramo que sale de "actual", calculamos su dirección y
        # con eso la posición del nodo al que llega.
        for t in salientes[actual]:
            if dir_in_actual is None:
                # "actual" es la raíz: no hay tramo anterior del que
                # partir, así que el primer tramo usa la dirección inicial
                # fija por convención, y su ref (plano XY para el Giro del
                # primer accesorio) se arma con _ref_inicial().
                dir_salida = np.array(direccion_inicial, dtype=float)
                dir_salida = dir_salida / np.linalg.norm(dir_salida)
                ref_salida = _ref_inicial(dir_salida)
            else:
                dir_salida, ref_salida = _direccion_salida(
                    tipo_actual, dir_in_actual, ref_in_actual, params_actual, t["params"]
                )
                dir_salida = dir_salida / np.linalg.norm(dir_salida)

            L = t["params"]["L"]
            destino = t["destino"]
            # Si el destino ya tiene posición, la red vuelve a un nodo ya
            # visitado (forma un lazo): sin este control el recorrido no
            # termina nunca.
            if destino in resultado:
                raise ValueError(
                    f"La red tiene un lazo: el tramo {actual} -> {destino} "
                    f"vuelve a un nodo ya visitado. La red debe ser un árbol."
                )
            # posición del siguiente nodo = posición actual + L en esa dirección
            pos_destino = pos_actual + dir_salida * L
            resultado[destino] = {"pos": pos_destino, "dir_in": dir_salida, "ref_in": ref_salida}
            pendientes.append(destino)  # lo encolamos para seguir desde ahí

    return resultado


if __name__ == "__main__":
    import sys
    import pprint

    sys.path.insert(0, ".")
    from mermaid_reader import leer_mermaid

    if len(sys.argv) < 2:
        print("Uso: python red.py archivo.md")
        sys.exit(1)

    nodos, tramos = leer_mermaid(sys.argv[1])
    posiciones = calcular_posiciones(nodos, tramos)
    for nid, info in posiciones.items():
        print(f"{nid} ({nodos[nid]['tipo']}): pos={info['pos']}, dir_in={info['dir_in']}")
