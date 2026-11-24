"""
mesh_preforma_web.py — interfaz gráfica de Mesh_preforma que corre en el NAVEGADOR.

Es el reemplazo de mesh_preforma_gui.py para equipos donde tkinter no se puede
cargar (por ejemplo, cuando el Control de aplicaciones de Windows bloquea
_tkinter). No necesita tkinter, Qt ni ninguna librería extra: usa solo la
biblioteca estándar de Python (servidor HTTP local) + numpy, y la ventana es
una página HTML (web_gui/index.html) que se abre en el navegador.

Hace lo mismo que la interfaz de escritorio:
  - editor del Mermaid con botones "+ Tramo recto / + Tee / + Codo / + Reducción";
  - diagrama 2D y traza 3D de los ejes, que se actualizan solos al escribir;
  - "Previsualizar malla" / "Mallar": genera la malla real y la muestra en 3D
    (se rota con el mouse; rueda = zoom; Mayús+arrastrar o botón derecho = mover);
  - "Guardar .malla": descarga el archivo .malla (con soldaduras);
  - abrir / guardar el Mermaid (.md).

Uso:
    uv run python mesh_preforma_web.py [archivo.md] [--puerto 8765] [--no-abrir]

El servidor solo escucha en 127.0.0.1 (esta computadora). Se detiene con Ctrl+C.
"""

import argparse
import contextlib
import io
import json
import os
import re
import sys
import tempfile
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(BASE, "src"))
sys.path.insert(0, os.path.join(BASE, "src", "parts"))

from mermaid_reader import (  # noqa: E402
    parsear_mermaid,
    PATRON_NODO_LIBRE,
    PATRON_NODO_PARAMS,
    PATRON_TRAMO,
    _parsear_params,
)
from red import calcular_posiciones  # noqa: E402
from ensamblador import ensamblar, fusionar_piezas  # noqa: E402
from exportador import escribir_malla  # noqa: E402
import tee  # noqa: E402

RUTA_HTML = os.path.join(BASE, "web_gui", "index.html")

PLANTILLA = """```mermaid
flowchart TD
1((Libre))
2["Tee:OD2=100,Alpha=90,H1=200,H2=200,H3=200,Giro=0"]
3["Codo:R=225,Angulo=90,Giro=90"]
4["Reduccion:OD1=150,OD2=100,L=250"]
5((Libre))
6((Libre))

1 -->|L=1000,OD=150| 2
2 -->|L=2000,OD=150,S| 3
3 -->|L=1000,OD=150| 4
4 -->|L=1000,OD=100| 5
2 -->|L=1000,OD=100,D| 6
```
"""

PLANTILLA_VACIA = """```mermaid
flowchart TD
1((Libre))
```
"""

LOCK = threading.Lock()
ULTIMA_MALLA = {"bytes": None, "nombre": "malla.malla"}
ARCHIVO_INICIAL = {"texto": None, "nombre": None}


# ------------------------------------------------------------ utilidades -----
def _fmt(valor):
    """Número sin decimales de más (1000 en vez de 1000.0)."""
    try:
        f = float(valor)
    except (TypeError, ValueError):
        return str(valor)
    return str(int(f)) if f.is_integer() else str(f)


def _etiqueta_params(params):
    partes = []
    for clave, valor in params.items():
        partes.append(clave if valor is True else f"{clave}={_fmt(valor)}")
    return ",".join(partes)


def _layout_mermaid(nodos, tramos):
    """Posición (nivel, fila) de cada nodo para dibujar el flowchart."""
    hijos = {n: [] for n in nodos}
    entrada = {n: 0 for n in nodos}
    for t in tramos:
        if t["origen"] not in hijos or t["destino"] not in hijos:
            continue
        hijos[t["origen"]].append(t["destino"])
        entrada[t["destino"]] = entrada.get(t["destino"], 0) + 1
    raices = [n for n in nodos if entrada.get(n, 0) == 0] or list(nodos)[:1]

    pos, visitados, contador_y = {}, set(), [0]

    def visitar(n, x):
        if n in visitados:
            return pos[n][1]
        visitados.add(n)
        desc = [h for h in hijos.get(n, []) if h not in visitados]
        if not desc:
            y = contador_y[0]
            contador_y[0] += 1
        else:
            ys = [visitar(h, x + 1) for h in desc]
            y = sum(ys) / len(ys)
        pos[n] = (x, y)
        return y

    for r in raices:
        visitar(r, 0)
    for n in nodos:
        if n not in pos:
            pos[n] = (0, contador_y[0])
            contador_y[0] += 1
    return pos


def _ordenar_anillo(pts):
    """Ordena por ángulo los puntos de un anillo casi plano (para dibujarlo
    como línea cerrada)."""
    centro = pts.mean(axis=0)
    rel = pts - centro
    _, _, vt = np.linalg.svd(rel, full_matrices=False)
    normal = vt[-1]
    e1 = rel[np.argmax(np.linalg.norm(rel, axis=1))]
    e1 = e1 - np.dot(e1, normal) * normal
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(normal, e1)
    return pts[np.argsort(np.arctan2(rel @ e2, rel @ e1))]


def _traza(nodos, tramos):
    """Esqueleto 3D de la red (posiciones de nodos y segmentos de tramo).
    Devuelve None si la red todavía no es válida."""
    try:
        posiciones = calcular_posiciones(nodos, tramos)
    except Exception as e:  # red incompleta mientras se escribe
        return None, str(e)
    pts = [{"id": nid, "tipo": nodos[nid]["tipo"], "pos": [float(v) for v in info["pos"]]}
           for nid, info in posiciones.items()]
    segs = [[[float(v) for v in posiciones[t["origen"]]["pos"]],
             [float(v) for v in posiciones[t["destino"]]["pos"]]]
            for t in tramos if t["origen"] in posiciones and t["destino"] in posiciones]
    return {"nodos": pts, "segmentos": segs}, None


# ---------------------------------------- operaciones sobre el texto Mermaid ---
def _lineas(texto):
    return texto.replace("\r\n", "\n").split("\n")


def ids_actuales(texto):
    ids = []
    for linea in _lineas(texto):
        m = PATRON_NODO_LIBRE.match(linea) or PATRON_NODO_PARAMS.match(linea)
        if m:
            ids.append(m.group(1))
    return ids


def siguiente_id_excluyendo(usados):
    n = 1
    while str(n) in usados:
        n += 1
    return str(n)


def od_sugerido(texto):
    od = 150
    for linea in _lineas(texto):
        m = PATRON_TRAMO.match(linea)
        if m:
            params = _parsear_params(m.group(2) or "")
            if "OD" in params and params["OD"] is not True:
                od = params["OD"]
    return od


def grado_nodo(texto, nodo_id):
    grado = 0
    for linea in _lineas(texto):
        m = PATRON_TRAMO.match(linea)
        if m:
            grado += (m.group(1) == nodo_id) + (m.group(3) == nodo_id)
    return grado


def od_entrante(texto, nodo_id):
    for linea in _lineas(texto):
        m = PATRON_TRAMO.match(linea)
        if m and m.group(3) == nodo_id:
            od = _parsear_params(m.group(2) or "").get("OD")
            if od is not None and od is not True:
                return od
    return None


def tramo_entrante(texto, nodo_id):
    for linea in _lineas(texto):
        m = PATRON_TRAMO.match(linea)
        if m and m.group(3) == nodo_id:
            return True
    return False


def insertar_lineas(texto, nuevas_decls, nuevas_tramos, reemplazo_nodo=None):
    """Igual que App._insertar_lineas de la GUI de escritorio: inserta las
    declaraciones tras el último nodo y los tramos tras el último tramo; con
    reemplazo_nodo convierte un extremo Libre ya conectado en el accesorio."""
    lineas = _lineas(texto)
    if reemplazo_nodo is not None:
        nodo_id, nueva_decl = reemplazo_nodo
        for i, linea in enumerate(lineas):
            m = PATRON_NODO_LIBRE.match(linea) or PATRON_NODO_PARAMS.match(linea)
            if m and m.group(1) == nodo_id:
                lineas[i] = nueva_decl
                break

    idx_ultimo_nodo = idx_ultimo_tramo = None
    idx_cierre = len(lineas)
    for i, linea in enumerate(lineas):
        if linea.strip().startswith("```") and i > 0:
            idx_cierre = i
            break
        if PATRON_NODO_LIBRE.match(linea) or PATRON_NODO_PARAMS.match(linea):
            idx_ultimo_nodo = i
        elif PATRON_TRAMO.match(linea):
            idx_ultimo_tramo = i

    punto_tramo = (idx_ultimo_tramo + 1) if idx_ultimo_tramo is not None else idx_cierre
    for linea in reversed(nuevas_tramos):
        lineas.insert(punto_tramo, linea)
    punto_nodo = (idx_ultimo_nodo + 1) if idx_ultimo_nodo is not None else 2
    for linea in reversed(nuevas_decls):
        lineas.insert(punto_nodo, linea)
    return "\n".join(lineas)


def resolver_origen(texto, origen):
    ids = ids_actuales(texto)
    if not ids:
        raise ValueError("No hay ningún nodo en el Mermaid todavía. Usá 'Nuevo archivo' o escribí al menos un nodo primero.")
    origen = (origen or "").strip()
    return origen if origen in ids else ids[-1]


MSJ_NODO_NO_DISPONIBLE = (
    "El nodo {o} ya tiene un tramo saliendo de él (no entrante), así que no es un "
    "extremo libre donde seguir agregando algo. Elegí en 'Agregar desde nodo' un extremo "
    "libre real (recién creado, sin tramos todavía, o con un único tramo entrante)."
)


def formulario(tipo, texto, origen):
    """Campos que hay que pedir para agregar `tipo` desde `origen`.
    Devuelve (titulo, campos) o lanza ValueError con el aviso a mostrar."""
    origen = resolver_origen(texto, origen)
    grado = grado_nodo(texto, origen)
    od_sug = od_sugerido(texto)

    if tipo == "tramo":
        if grado >= 1:
            raise ValueError(
                f"El nodo {origen} ya tiene un tramo. 'Tramo recto' solo se usa para iniciar la red "
                "desde un nodo todavía suelto. Para continuar desde acá agregá un accesorio "
                "(+ Codo / + Tee / + Reducción) — o si lo que querés es más longitud, editá el "
                "'L' del tramo existente directamente en el texto.")
        return "Agregar tramo recto", [
            ["L", "Longitud (L):", 1000], ["OD", "Diámetro (OD):", od_sug]]

    if grado >= 1 and not tramo_entrante(texto, origen):
        raise ValueError(MSJ_NODO_NO_DISPONIBLE.format(o=origen))

    entrada = [["L_in", "Longitud tramo de entrada:", 1000],
               ["OD_in", "Diámetro tramo de entrada:", od_sug]]
    if tipo == "tee":
        campos = [["OD2", "OD2 (derivación):", od_sug], ["Alpha", "Alpha (°):", 90],
                  ["H1", "H1:", 200], ["H2", "H2:", 200], ["H3", "H3:", 200],
                  ["Giro", "Giro (°):", 0],
                  ["L_recto", "Longitud tramo recto (S):", 1000],
                  ["L_deriv", "Longitud tramo derivación (D):", 1000]]
        return "Agregar Tee", (entrada + campos if grado == 0 else campos)
    if tipo == "codo":
        campos = [["R", "R (radio de curvatura):", 225], ["Angulo", "Ángulo del codo (°):", 90],
                  ["Giro", "Giro (°):", 0], ["L_out", "Longitud tramo de salida:", 1000]]
        return "Agregar Codo", (entrada + campos if grado == 0 else campos)
    if tipo == "reduccion":
        campos = [["OD1", "Diámetro de entrada (OD1):", od_sug],
                  ["L", "Longitud de la reducción:", 250],
                  ["OD2", "Diámetro de salida (OD2):", max(od_sug - 50, 50)],
                  ["L_out", "Longitud tramo de salida:", 1000]]
        return "Agregar Reducción", ([["L_in", "Longitud tramo de entrada:", 1000]] + campos
                                      if grado == 0 else campos)
    raise ValueError(f"Tipo desconocido: {tipo}")


def agregar(tipo, texto, origen, valores):
    """Agrega el elemento al texto. Devuelve (texto_nuevo, nuevo_origen, mensaje)."""
    titulo, _ = formulario(tipo, texto, origen)  # revalida el origen
    origen = resolver_origen(texto, origen)
    grado = grado_nodo(texto, origen)
    od_sug = od_sugerido(texto)
    try:
        v = {k: float(x) for k, x in valores.items()}
    except (TypeError, ValueError):
        raise ValueError("Todos los campos deben ser números.")

    usados = set(ids_actuales(texto))

    if tipo == "tramo":
        nuevo = siguiente_id_excluyendo(usados)
        t = insertar_lineas(texto, [f"{nuevo}((Libre))"],
                            [f"{origen} -->|L={_fmt(v['L'])},OD={_fmt(v['OD'])}| {nuevo}"])
        return t, nuevo, f"Tramo recto agregado: {origen} -> {nuevo}"

    if tipo == "tee":
        v["OD_recto"] = v["OD_in"] if grado == 0 else (od_entrante(texto, origen) or od_sug)
        if grado == 0:
            tee_id = siguiente_id_excluyendo(usados)
            usados.add(tee_id)
        else:
            tee_id = origen
        nodo_s = siguiente_id_excluyendo(usados)
        usados.add(nodo_s)
        nodo_d = siguiente_id_excluyendo(usados)
        decl_tee = (f'{tee_id}["Tee:OD2={_fmt(v["OD2"])},Alpha={_fmt(v["Alpha"])},'
                    f'H1={_fmt(v["H1"])},H2={_fmt(v["H2"])},H3={_fmt(v["H3"])},Giro={_fmt(v["Giro"])}"]')
        tramo_s = f'{tee_id} -->|L={_fmt(v["L_recto"])},OD={_fmt(v["OD_recto"])},S| {nodo_s}'
        tramo_d = f'{tee_id} -->|L={_fmt(v["L_deriv"])},OD={_fmt(v["OD2"])},D| {nodo_d}'
        decl_s, decl_d = f"{nodo_s}((Libre))", f"{nodo_d}((Libre))"
        if grado == 0:
            tramo_in = f'{origen} -->|L={_fmt(v["L_in"])},OD={_fmt(v["OD_in"])}| {tee_id}'
            t = insertar_lineas(texto, [decl_tee, decl_s, decl_d], [tramo_in, tramo_s, tramo_d])
        else:
            t = insertar_lineas(texto, [decl_s, decl_d], [tramo_s, tramo_d],
                                reemplazo_nodo=(tee_id, decl_tee))
        return t, nodo_s, (f"Tee en {tee_id}: recto -> {nodo_s} / derivación -> {nodo_d}. "
                           f"Origen actual: {nodo_s} (cambialo si querés seguir por la derivación {nodo_d}).")

    if tipo == "codo":
        v["OD_out"] = v["OD_in"] if grado == 0 else (od_entrante(texto, origen) or od_sug)
        if grado == 0:
            codo_id = siguiente_id_excluyendo(usados)
            usados.add(codo_id)
        else:
            codo_id = origen
        nodo_out = siguiente_id_excluyendo(usados)
        decl = f'{codo_id}["Codo:R={_fmt(v["R"])},Angulo={_fmt(v["Angulo"])},Giro={_fmt(v["Giro"])}"]'
        decl_out = f"{nodo_out}((Libre))"
        tramo_out = f'{codo_id} -->|L={_fmt(v["L_out"])},OD={_fmt(v["OD_out"])}| {nodo_out}'
        if grado == 0:
            tramo_in = f'{origen} -->|L={_fmt(v["L_in"])},OD={_fmt(v["OD_in"])}| {codo_id}'
            t = insertar_lineas(texto, [decl, decl_out], [tramo_in, tramo_out])
        else:
            t = insertar_lineas(texto, [decl_out], [tramo_out], reemplazo_nodo=(codo_id, decl))
        return t, nodo_out, f"Codo en {codo_id} -> {nodo_out}."

    if tipo == "reduccion":
        if grado == 0:
            red_id = siguiente_id_excluyendo(usados)
            usados.add(red_id)
        else:
            red_id = origen
        nodo_out = siguiente_id_excluyendo(usados)
        decl = f'{red_id}["Reduccion:OD1={_fmt(v["OD1"])},OD2={_fmt(v["OD2"])},L={_fmt(v["L"])}"]'
        decl_out = f"{nodo_out}((Libre))"
        tramo_out = f'{red_id} -->|L={_fmt(v["L_out"])},OD={_fmt(v["OD2"])}| {nodo_out}'
        if grado == 0:
            tramo_in = f'{origen} -->|L={_fmt(v["L_in"])},OD={_fmt(v["OD1"])}| {red_id}'
            t = insertar_lineas(texto, [decl, decl_out], [tramo_in, tramo_out])
        else:
            t = insertar_lineas(texto, [decl_out], [tramo_out], reemplazo_nodo=(red_id, decl))
        return t, nodo_out, f"Reducción en {red_id} -> {nodo_out}."

    raise ValueError(f"Tipo desconocido: {tipo}")


# ------------------------------------------------------------- endpoints ------
def api_info(d):
    """Diagrama 2D + traza 3D del texto actual (liviano, sin mallar)."""
    texto = d.get("texto", "")
    res = {"ids": ids_actuales(texto), "diagrama": None, "traza": None, "error": None}
    try:
        nodos, tramos = parsear_mermaid(texto)
    except Exception as e:
        res["error"] = str(e)
        return res
    if not nodos:
        return res
    lay = _layout_mermaid(nodos, tramos)
    res["diagrama"] = {
        "nodos": [{"id": n, "tipo": info["tipo"], "etiqueta": _etiqueta_params(info["params"]),
                   "x": lay[n][0], "y": lay[n][1]} for n, info in nodos.items()],
        "tramos": [{"origen": t["origen"], "destino": t["destino"],
                    "etiqueta": _etiqueta_params(t["params"])} for t in tramos],
    }
    res["traza"], res["error_traza"] = _traza(nodos, tramos)
    return res


def api_mallar(d):
    texto = d.get("texto", "")
    N = int(d.get("N", 24))
    if not 4 <= N <= 400:
        raise ValueError("N debe estar entre 4 y 400.")
    buf = io.StringIO()
    with LOCK:
        nodos, tramos = parsear_mermaid(texto)
        if not nodos:
            raise ValueError("No se encontró ningún nodo en el texto.")
        with contextlib.redirect_stdout(buf):
            piezas = ensamblar(nodos, tramos, N=N, verbose=True)
            posiciones = calcular_posiciones(nodos, tramos)
            nodos_f, elementos_f, zonas = fusionar_piezas(piezas)
            control = tee.control_malla(nodos_f, elementos_f)

        # .malla listo para descargar
        with tempfile.TemporaryDirectory() as tmp:
            ruta = os.path.join(tmp, "malla.malla")
            escribir_malla(ruta, nodos_f, elementos_f, mermaid_texto=texto, zonas_soldadura=zonas)
            with open(ruta, "rb") as f:
                ULTIMA_MALLA["bytes"] = f.read()
        base = os.path.splitext(ARCHIVO_INICIAL["nombre"] or "red")[0]
        ULTIMA_MALLA["nombre"] = base + ".malla"

    anillos = []
    for ids in zonas:
        anillo = _ordenar_anillo(nodos_f[np.asarray(ids, int)])
        centro = anillo.mean(axis=0)
        anillo = centro + (anillo - centro) * 1.08
        anillos.append(np.round(anillo, 3).ravel().tolist())

    traza, _ = _traza(nodos, tramos)
    return {
        "piezas": len(piezas), "N": N,
        "nodos": np.round(nodos_f, 3).ravel().tolist(),
        "elementos": np.asarray(elementos_f, int).ravel().tolist(),
        "anillos": anillos, "traza": traza,
        "control": {k: (float(v) if isinstance(v, (float, np.floating)) else int(v))
                    for k, v in control.items()},
        "n_zonas": len(zonas), "n_nodos_union": int(sum(len(z) for z in zonas)),
        "log": buf.getvalue(),
    }


def api_formulario(d):
    titulo, campos = formulario(d["tipo"], d.get("texto", ""), d.get("origen"))
    return {"titulo": titulo, "campos": campos}


def api_agregar(d):
    texto, nuevo, msj = agregar(d["tipo"], d.get("texto", ""), d.get("origen"), d.get("valores", {}))
    return {"texto": texto, "origen": nuevo, "mensaje": msj}


def api_inicio(_d):
    return {"texto": ARCHIVO_INICIAL["texto"] or PLANTILLA,
            "nombre": ARCHIVO_INICIAL["nombre"] or "red.md", "vacia": PLANTILLA_VACIA}


RUTAS_POST = {"/api/info": api_info, "/api/mallar": api_mallar, "/api/formulario": api_formulario,
              "/api/agregar": api_agregar, "/api/inicio": api_inicio}


class Manejador(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):  # sin ruido en consola
        pass

    def _enviar(self, codigo, cuerpo, tipo, extra=None):
        self.send_response(codigo)
        self.send_header("Content-Type", tipo)
        self.send_header("Content-Length", str(len(cuerpo)))
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(cuerpo)

    def do_GET(self):
        ruta = self.path.split("?")[0]
        if ruta in ("/", "/index.html"):
            with open(RUTA_HTML, "rb") as f:
                self._enviar(200, f.read(), "text/html; charset=utf-8")
        elif ruta == "/api/descargar":
            if ULTIMA_MALLA["bytes"] is None:
                self._enviar(404, "Todavía no se generó ninguna malla.".encode("utf-8"), "text/plain; charset=utf-8")
            else:
                self._enviar(200, ULTIMA_MALLA["bytes"], "application/octet-stream",
                             {"Content-Disposition": f'attachment; filename="{ULTIMA_MALLA["nombre"]}"'})
        elif ruta == "/favicon.ico":
            self._enviar(204, b"", "image/x-icon")
        else:
            self._enviar(404, b"no encontrado", "text/plain")

    def do_POST(self):
        fn = RUTAS_POST.get(self.path)
        if fn is None:
            self._enviar(404, b"no encontrado", "text/plain")
            return
        try:
            n = int(self.headers.get("Content-Length", 0))
            datos = json.loads(self.rfile.read(n) or b"{}")
            salida = {"ok": True, "datos": fn(datos)}
        except Exception as e:  # el mensaje llega al usuario tal cual
            salida = {"ok": False, "error": str(e) or e.__class__.__name__}
        self._enviar(200, json.dumps(salida).encode("utf-8"), "application/json; charset=utf-8")


def main():
    ap = argparse.ArgumentParser(description="Mesh_preforma — interfaz en el navegador")
    ap.add_argument("archivo", nargs="?", help="Mermaid (.md) para cargar al abrir")
    ap.add_argument("--puerto", type=int, default=8765)
    ap.add_argument("--no-abrir", action="store_true", help="no abrir el navegador solo")
    a = ap.parse_args()

    if a.archivo:
        with open(a.archivo, encoding="utf-8") as f:
            ARCHIVO_INICIAL["texto"] = f.read()
        ARCHIVO_INICIAL["nombre"] = os.path.basename(a.archivo)

    srv = ThreadingHTTPServer(("127.0.0.1", a.puerto), Manejador)
    url = f"http://127.0.0.1:{a.puerto}/"
    print(f"Mesh_preforma — interfaz web en {url}   (Ctrl+C para salir)")
    if not a.no_abrir:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nCerrado.")
    finally:
        srv.server_close()


if __name__ == "__main__":
    main()
