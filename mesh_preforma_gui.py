"""
mesh_preforma_gui.py — entorno gráfico de escritorio para armar la red de cañerías:

  - Editor de texto integrado para escribir/editar el Mermaid, con botones
    de ayuda ("Tramo recto", "Tee", "Codo") que agregan los nodos y tramos
    correspondientes automáticamente, pidiendo sólo los datos numéricos en
    un formulario.
  - Visualizador del diagrama Mermaid (2D) y traza 3D (ejes, liviana): se
    actualizan solos mientras escribís, sin apretar ningún botón.
  - "Previsualizar malla": genera la malla real (ensambla cada pieza y las
    funde) y la muestra en 3D, igual que "Mallar" pero sin el cartel de
    advertencia — para ver de entrada cómo va a quedar antes de guardar.
  - "Mallar": lo mismo, con el control de malla completo (avisa si hay
    problemas de topología u orientación antes de exportar).
  - "Guardar .malla...": exporta la malla ya generada a un archivo .malla
    (con el Mermaid de origen como encabezado).
  - "Nuevo archivo": limpia el editor para empezar otro caso desde cero.

Uso:
    uv run python mesh_preforma_gui.py [archivo.md]

    archivo.md : opcional, un Mermaid para cargar al abrir la ventana. Si
                 no se pasa nada, arranca con una plantilla de ejemplo.
"""

import contextlib
import io
import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

import numpy as np
import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, "src")
sys.path.insert(0, "src/parts")

from mermaid_reader import (
    parsear_mermaid,
    PATRON_NODO_LIBRE,
    PATRON_NODO_PARAMS,
    PATRON_TRAMO,
    _parsear_params,
)
from red import calcular_posiciones
from ensamblador import ensamblar, fusionar_piezas
from exportador import escribir_malla
import graficador
import tee


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

COLORES_TIPO = {
    "Libre": "#d9d9d9",
    "Tubo": "#cfe8ff",
    "Codo": "#ffd9a0",
    "Reduccion": "#c9f2c9",
    "Tee": "#f2c9e8",
}


def _fmt(valor):
    """Formatea un número sin decimales de más (1000 en vez de 1000.0)."""
    try:
        f = float(valor)
    except (TypeError, ValueError):
        return str(valor)
    if f.is_integer():
        return str(int(f))
    return str(f)


def _etiqueta_params(params):
    """'OD2=100,Alpha=90,S' a partir de un dict de params (flags True se
    muestran solo por su nombre, sin '=True')."""
    partes = []
    for clave, valor in params.items():
        if valor is True:
            partes.append(clave)
        else:
            partes.append(f"{clave}={_fmt(valor)}")
    return ",".join(partes)


def _layout_mermaid(nodos, tramos):
    """Ubica cada nodo en una grilla (nivel = profundidad desde una raíz,
    posición vertical = orden tipo árbol) para poder dibujar el flowchart.
    Devuelve {nodo_id: (x, y)}."""
    hijos = {n: [] for n in nodos}
    entrada = {n: 0 for n in nodos}
    for t in tramos:
        if t["origen"] not in hijos or t["destino"] not in hijos:
            continue
        hijos[t["origen"]].append(t["destino"])
        entrada[t["destino"]] = entrada.get(t["destino"], 0) + 1

    raices = [n for n in nodos if entrada.get(n, 0) == 0]
    if not raices:
        raices = list(nodos)[:1]

    pos = {}
    visitados = set()
    contador_y = [0]

    def visitar(n, x):
        if n in visitados:
            return pos[n][1]
        visitados.add(n)
        descendientes = [h for h in hijos.get(n, []) if h not in visitados]
        if not descendientes:
            y = contador_y[0]
            contador_y[0] += 1
        else:
            ys = [visitar(h, x + 1) for h in descendientes]
            y = sum(ys) / len(ys)
        pos[n] = (x, y)
        return y

    for r in raices:
        visitar(r, 0)
    for n in nodos:
        if n not in pos:
            y = contador_y[0]
            contador_y[0] += 1
            pos[n] = (0, y)
    return pos


class App(tk.Tk):
    def __init__(self, archivo_inicial=None):
        super().__init__()
        self.title("Mesh Preforma — editor de red de cañerías")
        self.geometry("1360x860")

        # estado del último mallado hecho (None hasta apretar "Mallar"; se
        # invalida cada vez que se edita o se vuelve a previsualizar, para
        # no permitir guardar un .malla que ya no corresponde al texto
        # actual del editor)
        self.piezas = None
        self.nodos_f = None
        self.elementos_f = None
        self.zonas_soldadura = None
        self.archivo_actual = None  # ruta del .md cargado/guardado (o None)
        self._despues_auto_refresco = None  # id de after() pendiente (debounce)

        self._armar_ui()

        # refresco automático (diagrama 2D + traza 3D liviana) mientras se
        # escribe, sin tener que apretar ningún botón
        self.editor.edit_modified(False)
        self.editor.bind("<<Modified>>", self._on_editor_modificado)

        if archivo_inicial and os.path.exists(archivo_inicial):
            self._cargar_archivo(archivo_inicial)
        else:
            self.editor.insert("1.0", PLANTILLA)
            self.editor.edit_modified(False)
            self._actualizar_combo_origen()
            self._dibujar_diagrama()
            self._dibujar_traza_automatica()

    # ------------------------------------------------------------ UI -----
    def _armar_ui(self):
        barra = ttk.Frame(self)
        barra.pack(side="top", fill="x", padx=6, pady=(4, 0))

        ttk.Button(barra, text="Nuevo archivo", command=self._nuevo_archivo).pack(side="left", padx=2)
        ttk.Button(barra, text="Abrir...", command=self._abrir).pack(side="left", padx=2)
        ttk.Button(barra, text="Guardar Mermaid...", command=self._guardar_mermaid).pack(side="left", padx=2)
        ttk.Separator(barra, orient="vertical").pack(side="left", fill="y", padx=6)

        ttk.Button(barra, text="↶ Deshacer (Ctrl+Z)", command=self._deshacer).pack(side="left", padx=2)
        ttk.Button(barra, text="↷ Rehacer (Ctrl+Y)", command=self._rehacer).pack(side="left", padx=2)
        ttk.Separator(barra, orient="vertical").pack(side="left", fill="y", padx=6)

        ttk.Label(barra, text="N (divisiones de perímetro):").pack(side="left")
        self.var_n = tk.IntVar(value=24)
        ttk.Spinbox(barra, from_=4, to=400, increment=2, textvariable=self.var_n, width=6).pack(
            side="left", padx=4
        )

        ttk.Button(barra, text="Previsualizar malla", command=self._previsualizar).pack(
            side="left", padx=(16, 2)
        )
        ttk.Button(barra, text="Mallar", command=self._mallar).pack(side="left", padx=2)
        ttk.Button(barra, text="Guardar .malla...", command=self._guardar_malla).pack(
            side="left", padx=(16, 2)
        )

        # --- segunda barra: ayudas de escritura del mermaid ---
        barra2 = ttk.Frame(self)
        barra2.pack(side="top", fill="x", padx=6, pady=(4, 4))

        ttk.Label(barra2, text="Agregar desde nodo:").pack(side="left")
        self.var_origen = tk.StringVar()
        self.combo_origen = ttk.Combobox(barra2, textvariable=self.var_origen, width=6, state="readonly")
        self.combo_origen.pack(side="left", padx=4)

        ttk.Button(barra2, text="+ Tramo recto", command=self._agregar_tramo_recto).pack(
            side="left", padx=(12, 2)
        )
        ttk.Button(barra2, text="+ Tee", command=self._agregar_tee).pack(side="left", padx=2)
        ttk.Button(barra2, text="+ Codo", command=self._agregar_codo).pack(side="left", padx=2)
        ttk.Button(barra2, text="+ Reducción", command=self._agregar_reduccion).pack(side="left", padx=2)

        panel = ttk.PanedWindow(self, orient="horizontal")
        panel.pack(side="top", fill="both", expand=True)

        # --- panel izquierdo: editor de texto del Mermaid ---
        marco_editor = ttk.Frame(panel)
        ttk.Label(marco_editor, text="Mermaid (editable):").pack(anchor="w")
        self.editor = scrolledtext.ScrolledText(
            marco_editor, wrap="none", width=55, font=("Consolas", 10),
            undo=True, autoseparators=True, maxundo=-1,
        )
        self.editor.pack(fill="both", expand=True)
        self.editor.bind("<Control-z>", self._deshacer)
        self.editor.bind("<Control-Z>", self._deshacer)
        self.editor.bind("<Control-y>", self._rehacer)
        self.editor.bind("<Control-Y>", self._rehacer)
        self.editor.bind("<Control-Shift-Z>", self._rehacer)
        self.editor.bind("<Control-Shift-z>", self._rehacer)
        panel.add(marco_editor, weight=1)

        # --- panel derecho: pestañas (diagrama 2D / vista 3D) + log -----
        marco_vista = ttk.Frame(panel)
        self.notebook = ttk.Notebook(marco_vista)
        self.notebook.pack(fill="both", expand=True)

        # pestaña: diagrama del mermaid (2D)
        tab_diagrama = ttk.Frame(self.notebook)
        barra_diagrama = ttk.Frame(tab_diagrama)
        barra_diagrama.pack(side="top", fill="x")
        ttk.Button(barra_diagrama, text="Actualizar diagrama", command=self._dibujar_diagrama).pack(
            side="left", padx=4, pady=2
        )
        self.fig_diagrama = Figure(figsize=(7, 6))
        self.ax_diagrama = self.fig_diagrama.add_subplot(111)
        self.canvas_diagrama = FigureCanvasTkAgg(self.fig_diagrama, master=tab_diagrama)
        self.canvas_diagrama.get_tk_widget().pack(fill="both", expand=True)
        self.notebook.add(tab_diagrama, text="Diagrama Mermaid")

        # pestaña: vista 3D (traza / malla)
        tab_3d = ttk.Frame(self.notebook)
        self.fig = Figure(figsize=(7, 6))
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.canvas = FigureCanvasTkAgg(self.fig, master=tab_3d)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        barra_tk = NavigationToolbar2Tk(self.canvas, tab_3d)
        barra_tk.update()
        self.notebook.add(tab_3d, text="Vista 3D")

        self.notebook.bind("<<NotebookTabChanged>>", lambda e: self._dibujar_diagrama())

        ttk.Label(marco_vista, text="Estado:").pack(anchor="w")
        self.log = scrolledtext.ScrolledText(marco_vista, height=9, font=("Consolas", 9), state="disabled")
        self.log.pack(fill="x")

        panel.add(marco_vista, weight=2)

    # --------------------------------------------------------- deshacer -----
    def _reemplazar_texto_editor(self, nuevo_texto, reiniciar_historial=False):
        """Reemplaza todo el contenido del editor. El borrado + la escritura
        quedan agrupados como un solo paso de deshacer (un Ctrl+Z los revierte
        juntos), en vez de dejar al usuario con el editor vacío a mitad de
        camino."""
        self.editor.edit_separator()
        self.editor.delete("1.0", "end")
        self.editor.insert("1.0", nuevo_texto)
        self.editor.edit_separator()
        if reiniciar_historial:
            # al cargar un archivo nuevo o empezar un caso desde cero no
            # tiene sentido poder "deshacer" hasta el archivo anterior
            self.editor.edit_reset()

    def _deshacer(self, event=None):
        try:
            self.editor.edit_undo()
        except tk.TclError:
            self._escribir_log("No hay nada para deshacer.")
        else:
            self._actualizar_combo_origen()
            self._dibujar_diagrama()
            self._invalidar_mallado()
        return "break"

    def _rehacer(self, event=None):
        try:
            self.editor.edit_redo()
        except tk.TclError:
            self._escribir_log("No hay nada para rehacer.")
        else:
            self._actualizar_combo_origen()
            self._dibujar_diagrama()
            self._invalidar_mallado()
        return "break"

    def _escribir_log(self, texto):
        if not texto:
            return
        self.log.configure(state="normal")
        self.log.insert("end", texto.rstrip("\n") + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def _limpiar_log(self):
        self.log.configure(state="normal")
        self.log.delete("1.0", "end")
        self.log.configure(state="disabled")

    def _ejecutar_con_log(self, fn, *args, **kwargs):
        """Corre fn capturando lo que imprima por stdout (los avisos de
        mermaid_reader/ensamblador) y lo manda al log de la ventana en vez
        de perderlo en la consola."""
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            resultado = fn(*args, **kwargs)
        self._escribir_log(buffer.getvalue())
        return resultado

    # -------------------------------------------------------- archivo -----
    def _abrir(self):
        ruta = filedialog.askopenfilename(
            filetypes=[("Mermaid", "*.md *.mmd"), ("Todos los archivos", "*.*")]
        )
        if ruta:
            self._cargar_archivo(ruta)

    def _cargar_archivo(self, ruta):
        with open(ruta, "r", encoding="utf-8") as f:
            contenido = f.read()
        self._reemplazar_texto_editor(contenido, reiniciar_historial=True)
        self.archivo_actual = ruta
        self.title(f"Mesh Preforma — {os.path.basename(ruta)}")
        self._invalidar_mallado()
        self._actualizar_combo_origen()
        self._dibujar_diagrama()
        self._escribir_log(f"Cargado: {ruta}")

    def _guardar_mermaid(self):
        inicial = os.path.basename(self.archivo_actual) if self.archivo_actual else "red.md"
        ruta = filedialog.asksaveasfilename(
            defaultextension=".md", filetypes=[("Mermaid", "*.md")], initialfile=inicial
        )
        if not ruta:
            return
        with open(ruta, "w", encoding="utf-8") as f:
            f.write(self.editor.get("1.0", "end-1c"))
        self.archivo_actual = ruta
        self.title(f"Mesh Preforma — {os.path.basename(ruta)}")
        self._escribir_log(f"Mermaid guardado: {ruta}")

    def _nuevo_archivo(self):
        if not messagebox.askyesno(
            "Nuevo archivo", "¿Descartar el contenido actual del editor y empezar un caso nuevo?"
        ):
            return
        self._reemplazar_texto_editor(PLANTILLA_VACIA, reiniciar_historial=True)
        self.archivo_actual = None
        self.title("Mesh Preforma — nuevo caso")
        self._invalidar_mallado()
        self._actualizar_combo_origen()
        self._dibujar_diagrama()
        self._escribir_log("Nuevo caso iniciado. Usá los botones '+ Tramo recto / Tee / Codo' para armarlo, o escribí directo en el editor.")

    # --------------------------------------------------- previsualizar -----
    def _leer_editor(self):
        """Parsea el contenido actual del editor. Devuelve (nodos, tramos) o
        None si hay un error (ya mostrado al usuario con un messagebox)."""
        texto = self.editor.get("1.0", "end-1c")
        try:
            nodos, tramos = self._ejecutar_con_log(parsear_mermaid, texto)
        except Exception as e:
            messagebox.showerror("Error al leer el Mermaid", str(e))
            self._escribir_log(f"ERROR al leer el Mermaid: {e}")
            return None
        if not nodos:
            messagebox.showwarning("Mermaid vacío", "No se encontró ningún nodo en el texto.")
            return None
        return nodos, tramos

    def _invalidar_mallado(self):
        """El mallado que haya quedado de una corrida anterior ya no
        corresponde al texto actual del editor (se cambió o se recargó algo
        nuevo): hay que volver a mallar antes de poder guardar un .malla."""
        self.piezas = None
        self.nodos_f = self.elementos_f = self.zonas_soldadura = None

    def _generar_malla(self, titulo, avisar_problemas):
        """Hace el trabajo pesado de verdad (ensamblar + fusionar) y dibuja
        la malla real (tubos, no solo ejes) en la pestaña 3D. La usan tanto
        'Previsualizar' como 'Mallar' — así lo que ves al previsualizar es
        exactamente la malla que se va a generar, no un esqueleto aparte."""
        resultado = self._leer_editor()
        if resultado is None:
            return False
        nodos, tramos = resultado
        N = self.var_n.get()
        self._limpiar_log()
        try:
            piezas = self._ejecutar_con_log(ensamblar, nodos, tramos, N=N, verbose=True)
            posiciones = calcular_posiciones(nodos, tramos)
            nodos_f, elementos_f, zonas = fusionar_piezas(piezas)
            control = tee.control_malla(nodos_f, elementos_f)
        except Exception as e:
            messagebox.showerror("Error al mallar", str(e))
            self._escribir_log(f"ERROR al mallar: {e}")
            return False

        self.piezas = piezas
        self.nodos_f, self.elementos_f, self.zonas_soldadura = nodos_f, elementos_f, zonas

        self.ax.clear()
        graficador.graficar_malla(nodos_f, elementos_f, ax=self.ax, mostrar=False, colores={0: "#ffcf80"})
        graficador.graficar_soldaduras(nodos_f, zonas, ax=self.ax, mostrar=False)
        graficador.graficar_ejes(posiciones, tramos, nodos, ax=self.ax, mostrar=False)
        self._ajustar_zoom(nodos_f)
        self.ax.set_title(f"{titulo} — {len(piezas)} piezas, N={N}")
        self.canvas.draw()
        self.notebook.select(1)

        self._escribir_log(f"{titulo}: {len(piezas)} piezas, N={N}")
        self._escribir_log(f"Control de malla: {control}")
        self._escribir_log(
            f"Zonas de soldadura: {len(zonas)} ({sum(len(z) for z in zonas)} nodos de union)"
        )
        if avisar_problemas and (control["mas_de_2"] or control["orientacion_inconsistente"]):
            messagebox.showwarning(
                "Revisar malla",
                "El control de malla encontró problemas (ver el panel de Estado): "
                f"mas_de_2={control['mas_de_2']}, "
                f"orientacion_inconsistente={control['orientacion_inconsistente']}.",
            )
        return True

    def _previsualizar(self):
        # Antes solo mostraba el esqueleto (ejes); ahora genera la malla real
        # y muestra eso, para ver de entrada el resultado tal como va a
        # quedar — sin el cartel de advertencia (que sí aparece al mallar),
        # para no interrumpir mientras todavía estás iterando.
        self._generar_malla("Previsualización de la malla", avisar_problemas=False)

    # ----------------------------------------------------------- mallar -----
    def _mallar(self):
        self._generar_malla("Malla fusionada", avisar_problemas=True)

    def _ajustar_zoom(self, puntos):
        puntos = np.asarray(puntos, dtype=float).reshape(-1, 3)
        lo, hi = puntos.min(axis=0), puntos.max(axis=0)
        centro = (lo + hi) / 2
        radio = max((hi - lo).max() / 2, 1e-9) * 1.1
        self.ax.set_xlim(centro[0] - radio, centro[0] + radio)
        self.ax.set_ylim(centro[1] - radio, centro[1] + radio)
        self.ax.set_zlim(centro[2] - radio, centro[2] + radio)
        self.ax.set_box_aspect((1, 1, 1))

    # -------------------------------------------------- refresco automático -----
    def _on_editor_modificado(self, event=None):
        """Se dispara con cada cambio en el editor (tipeo, undo/redo, los
        botones de ayuda, cargar un archivo...). Tkinter solo manda este
        evento una vez por cambio real si se resetea la bandera; por eso
        se resetea acá mismo antes de programar el refresco."""
        if not self.editor.edit_modified():
            return
        self.editor.edit_modified(False)
        self._invalidar_mallado()
        self._programar_auto_refresco()

    def _programar_auto_refresco(self):
        """Debounce: si hay un refresco ya programado lo cancela y arranca
        de nuevo, para no recalcular en cada tecla — recién cuando pasan
        ~500ms sin cambios."""
        if self._despues_auto_refresco is not None:
            self.after_cancel(self._despues_auto_refresco)
        self._despues_auto_refresco = self.after(500, self._auto_refrescar)

    def _auto_refrescar(self):
        self._despues_auto_refresco = None
        self._actualizar_combo_origen()
        self._dibujar_diagrama()
        self._dibujar_traza_automatica()

    def _dibujar_traza_automatica(self):
        """Versión liviana (sin mallar) para el refresco automático: solo
        el esqueleto de ejes, igual que antes mostraba 'Previsualizar'. No
        interrumpe al usuario con mensajes de error mientras escribe —
        eso lo sigue haciendo 'Previsualizar malla' cuando lo aprieta a
        propósito. Tampoco toca self.piezas/nodos_f: el refresco automático
        nunca invalida ni reemplaza una malla ya generada a propósito."""
        texto = self.editor.get("1.0", "end-1c")
        try:
            nodos, tramos = parsear_mermaid(texto)
            if not nodos:
                return
            posiciones = calcular_posiciones(nodos, tramos)
        except Exception:
            return
        self.ax.clear()
        graficador.graficar_ejes(posiciones, tramos, nodos, ax=self.ax, mostrar=False)
        self._ajustar_zoom([info["pos"] for info in posiciones.values()])
        self.ax.set_title("Traza (ejes) — vista automática, sin mallar")
        self.canvas.draw()

    # ----------------------------------------------------- guardar .malla -----
    def _guardar_malla(self):
        if self.nodos_f is None:
            messagebox.showwarning(
                "Todavía no mallaste", "Primero hacé click en 'Mallar' para generar la malla."
            )
            return
        inicial = "red.malla"
        if self.archivo_actual:
            inicial = os.path.splitext(os.path.basename(self.archivo_actual))[0] + ".malla"
        ruta = filedialog.asksaveasfilename(
            defaultextension=".malla", filetypes=[("Malla", "*.malla")], initialfile=inicial
        )
        if not ruta:
            return
        mermaid_texto = self.editor.get("1.0", "end-1c")
        escribir_malla(
            ruta,
            self.nodos_f,
            self.elementos_f,
            mermaid_texto=mermaid_texto,
            zonas_soldadura=self.zonas_soldadura,
        )
        self._escribir_log(f"Malla exportada: {ruta}")
        messagebox.showinfo("Guardado", f"Malla exportada a:\n{ruta}")

    # ------------------------------------------------- diagrama mermaid -----
    def _dibujar_diagrama(self):
        """Dibuja el flowchart (cajas + flechas) tal como está el editor
        ahora mismo. Si el texto tiene un error de sintaxis no interrumpe
        al usuario con un mensaje — eso ya lo hacen 'Previsualizar'/'Mallar'
        — simplemente no actualiza el dibujo."""
        texto = self.editor.get("1.0", "end-1c")
        try:
            nodos, tramos = parsear_mermaid(texto)
        except Exception:
            return

        self.ax_diagrama.clear()
        self.ax_diagrama.axis("off")
        if not nodos:
            self.canvas_diagrama.draw()
            return

        # profundidad (x del layout) -> eje vertical, hacia abajo, como un
        # "flowchart TD" de verdad; orden entre hermanos (y del layout) ->
        # eje horizontal, para separar las ramas.
        pos = _layout_mermaid(nodos, tramos)
        ancho_caja, alto_caja = 1.8, 0.7
        dx, dy = 2.8, 1.7
        coords = {nid: (y * dx, -x * dy) for nid, (x, y) in pos.items()}

        for t in tramos:
            if t["origen"] not in coords or t["destino"] not in coords:
                continue
            x0, y0 = coords[t["origen"]]
            x1, y1 = coords[t["destino"]]
            self.ax_diagrama.annotate(
                "",
                xy=(x1, y1 + alto_caja / 2),
                xytext=(x0, y0 - alto_caja / 2),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.3, shrinkA=2, shrinkB=2),
            )
            etiqueta = _etiqueta_params(t["params"])
            if etiqueta:
                xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
                self.ax_diagrama.text(
                    xm + 0.12, ym, etiqueta, fontsize=7, color="#555555", ha="left", va="center"
                )

        for nid, (cx, cy) in coords.items():
            tipo = nodos[nid]["tipo"]
            color = COLORES_TIPO.get(tipo, "#e8e8e8")
            caja = FancyBboxPatch(
                (cx - ancho_caja / 2, cy - alto_caja / 2),
                ancho_caja,
                alto_caja,
                boxstyle="round,pad=0.03",
                linewidth=1.0,
                edgecolor="#333333",
                facecolor=color,
            )
            self.ax_diagrama.add_patch(caja)
            resumen_params = _etiqueta_params(nodos[nid]["params"])
            texto_nodo = f"{nid}: {tipo}"
            if resumen_params:
                texto_nodo += f"\n{resumen_params}"
            self.ax_diagrama.text(cx, cy, texto_nodo, fontsize=7, ha="center", va="center")

        xs = [c[0] for c in coords.values()]
        ys = [c[1] for c in coords.values()]
        self.ax_diagrama.set_xlim(min(xs) - ancho_caja, max(xs) + ancho_caja * 1.8)
        self.ax_diagrama.set_ylim(min(ys) - alto_caja, max(ys) + alto_caja)
        self.ax_diagrama.set_aspect("auto")
        self.canvas_diagrama.draw()

    # ------------------------------------------------ ayudas de escritura -----
    def _ids_actuales(self):
        """ids de los nodos declarados en el editor, en el orden en que
        aparecen (sin parsear tramos — alcanza con las líneas de nodo)."""
        texto = self.editor.get("1.0", "end-1c")
        ids = []
        for linea in texto.splitlines():
            m = PATRON_NODO_LIBRE.match(linea) or PATRON_NODO_PARAMS.match(linea)
            if m:
                ids.append(m.group(1))
        return ids

    def _siguiente_id_excluyendo(self, usados):
        n = 1
        while str(n) in usados:
            n += 1
        return str(n)

    def _siguiente_id(self):
        return self._siguiente_id_excluyendo(set(self._ids_actuales()))

    def _actualizar_combo_origen(self, seleccionar=None):
        ids = self._ids_actuales()
        self.combo_origen["values"] = ids
        if seleccionar and seleccionar in ids:
            self.var_origen.set(seleccionar)
        elif self.var_origen.get() not in ids and ids:
            self.var_origen.set(ids[-1])
        elif not ids:
            self.var_origen.set("")

    def _fijar_nodo_origen(self, nodo):
        self._actualizar_combo_origen(seleccionar=nodo)

    def _nodo_origen_actual(self):
        self._actualizar_combo_origen()
        ids = self._ids_actuales()
        if not ids:
            return None
        seleccion = self.var_origen.get().strip()
        if seleccion and seleccion in ids:
            return seleccion
        return ids[-1]

    def _od_sugerido(self):
        """Toma el OD del último tramo escrito, para no tener que
        retipearlo cada vez (si no hay ninguno todavía, 150 por defecto)."""
        texto = self.editor.get("1.0", "end-1c")
        od = 150
        for linea in texto.splitlines():
            m = PATRON_TRAMO.match(linea)
            if m:
                params = _parsear_params(m.group(2))
                if "OD" in params and params["OD"] is not True:
                    od = params["OD"]
        return od

    def _grado_nodo(self, nodo_id):
        """Cantidad de tramos que tocan a nodo_id (como origen o como
        destino, sumado). El modelo exige un grado fijo según el tipo:
        Libre=1, Codo=2, Reduccion=2, Tee=3 (ver GRADOS_ESPERADOS en
        red.py) — por eso las ayudas de escritura necesitan saber el grado
        actual antes de decidir si hay que crear un nodo nuevo o convertir
        el existente."""
        texto = self.editor.get("1.0", "end-1c")
        grado = 0
        for linea in texto.splitlines():
            m = PATRON_TRAMO.match(linea)
            if m:
                if m.group(1) == nodo_id:
                    grado += 1
                if m.group(3) == nodo_id:
                    grado += 1
        return grado

    def _tramo_entrante(self, nodo_id):
        """True si nodo_id es destino de algún tramo (tiene algo entrando).
        Un extremo libre "normal" (el que va quedando al final de cada
        pieza agregada) siempre cae en este caso cuando ya tiene grado >= 1;
        si en cambio el grado viene de ser origen de un tramo (alguien
        seleccionó a mano un nodo que ya tiene una salida), no es un punto
        válido para seguir agregando algo ahí."""
        texto = self.editor.get("1.0", "end-1c")
        for linea in texto.splitlines():
            m = PATRON_TRAMO.match(linea)
            if m and m.group(3) == nodo_id:
                return True
        return False

    def _pedir_valores(self, titulo, campos):
        """Ventana modal simple con un Entry por campo. `campos` es una
        lista de (clave, etiqueta, valor_por_defecto). Devuelve un dict
        {clave: texto} o None si se cancela."""
        ventana = tk.Toplevel(self)
        ventana.title(titulo)
        ventana.transient(self)
        ventana.grab_set()
        ventana.resizable(False, False)

        entradas = {}
        for i, (clave, etiqueta, defecto) in enumerate(campos):
            ttk.Label(ventana, text=etiqueta).grid(row=i, column=0, sticky="w", padx=8, pady=3)
            var = tk.StringVar(value=_fmt(defecto))
            entrada = ttk.Entry(ventana, textvariable=var, width=14)
            entrada.grid(row=i, column=1, padx=8, pady=3)
            entradas[clave] = var

        estado = {"ok": False}

        def _ok(event=None):
            estado["ok"] = True
            ventana.destroy()

        def _cancelar(event=None):
            ventana.destroy()

        fila_botones = len(campos)
        marco_botones = ttk.Frame(ventana)
        marco_botones.grid(row=fila_botones, column=0, columnspan=2, pady=(6, 8))
        ttk.Button(marco_botones, text="Aceptar", command=_ok).pack(side="left", padx=4)
        ttk.Button(marco_botones, text="Cancelar", command=_cancelar).pack(side="left", padx=4)

        ventana.bind("<Return>", _ok)
        ventana.bind("<Escape>", _cancelar)
        list(entradas.values())[0] and ventana.focus_force()

        ventana.wait_window()
        if not estado["ok"]:
            return None
        return {clave: var.get().strip() for clave, var in entradas.items()}

    def _insertar_lineas(self, nuevas_decls, nuevas_tramos, reemplazo_nodo=None):
        """Inserta declaraciones de nodo después de la última línea de nodo
        existente, y tramos después del último tramo existente (o antes del
        cierre ``` si todavía no hay ninguno de alguno de los dos).

        Si se pasa reemplazo_nodo=(nodo_id, nueva_declaracion), además
        reemplaza en el lugar la línea de declaración de ese nodo. Se usa
        para "convertir" un extremo libre ya conectado en el accesorio
        nuevo (Codo/Tee/Reducción), en vez de dejarlo colgado como "Libre"
        con 2 o más tramos (que el chequeo de topología de red.py rechaza:
        un nodo Libre tiene que tener exactamente 1 tramo)."""
        texto = self.editor.get("1.0", "end-1c")
        lineas = texto.split("\n")

        if reemplazo_nodo is not None:
            nodo_id, nueva_decl = reemplazo_nodo
            for i, linea in enumerate(lineas):
                m = PATRON_NODO_LIBRE.match(linea) or PATRON_NODO_PARAMS.match(linea)
                if m and m.group(1) == nodo_id:
                    lineas[i] = nueva_decl
                    break

        idx_ultimo_nodo = None
        idx_ultimo_tramo = None
        idx_cierre = len(lineas)
        for i, linea in enumerate(lineas):
            s = linea.strip()
            if s.startswith("```") and i > 0:
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

        self._reemplazar_texto_editor("\n".join(lineas))
        self._invalidar_mallado()

    def _origen_o_aviso(self):
        """Resuelve el nodo de origen actual, o avisa y devuelve None si
        todavía no hay ningún nodo."""
        origen = self._nodo_origen_actual()
        if origen is None:
            messagebox.showwarning(
                "Sin nodo de origen",
                "No hay ningún nodo en el Mermaid todavía. Usá 'Nuevo archivo' o escribí al menos un nodo primero.",
            )
        return origen

    def _agregar_tramo_recto(self):
        # El tramo recto es solo para INICIAR la red (el primer tramo que
        # sale de un nodo todavía suelto, grado 0). Para seguir desde un
        # extremo ya conectado hay que agregar un accesorio (Codo/Tee/
        # Reducción): un tramo recto atrás de otro no tiene sentido como
        # pieza aparte (sería solamente más longitud del mismo tramo), así
        # que no se ofrece acá para no dejar un nodo "Libre" con 2 tramos.
        origen = self._origen_o_aviso()
        if origen is None:
            return
        if self._grado_nodo(origen) >= 1:
            messagebox.showwarning(
                "No corresponde acá",
                f"El nodo {origen} ya tiene un tramo. 'Tramo recto' solo se usa para iniciar la red "
                "desde un nodo todavía suelto. Para continuar desde acá agregá un accesorio "
                "(+ Codo / + Tee / + Reducción) — o si lo que querés es más longitud, editá el "
                "'L' del tramo existente directamente en el texto.",
            )
            return

        datos = self._pedir_valores(
            "Agregar tramo recto",
            [("L", "Longitud (L):", 1000), ("OD", "Diámetro (OD):", self._od_sugerido())],
        )
        if datos is None:
            return
        try:
            L, OD = float(datos["L"]), float(datos["OD"])
        except ValueError:
            messagebox.showerror("Valor inválido", "Longitud y diámetro deben ser números.")
            return

        nuevo = self._siguiente_id()
        self._insertar_lineas(
            [f"{nuevo}((Libre))"], [f"{origen} -->|L={_fmt(L)},OD={_fmt(OD)}| {nuevo}"]
        )
        self._fijar_nodo_origen(nuevo)
        self._dibujar_diagrama()
        self._escribir_log(f"Tramo recto agregado: {origen} -> {nuevo} (L={_fmt(L)}, OD={_fmt(OD)})")

    def _advertir_nodo_no_disponible(self, origen):
        messagebox.showwarning(
            "Nodo no disponible",
            f"El nodo {origen} ya tiene un tramo saliendo de él (no entrante), así que no es un "
            "extremo libre donde seguir agregando algo. Elegí en el combo 'Agregar desde nodo' un "
            "extremo libre real (recién creado, sin tramos todavía, o con un único tramo entrante).",
        )

    def _agregar_tee(self):
        origen = self._origen_o_aviso()
        if origen is None:
            return
        grado = self._grado_nodo(origen)
        if grado >= 1 and not self._tramo_entrante(origen):
            self._advertir_nodo_no_disponible(origen)
            return

        od_sug = self._od_sugerido()
        campos_tee = [
            ("OD2", "OD2 (derivación):", od_sug),
            ("Alpha", "Alpha (°):", 90),
            ("H1", "H1:", 200),
            ("H2", "H2:", 200),
            ("H3", "H3:", 200),
            ("Giro", "Giro (°):", 0),
            ("L_recto", "Longitud tramo recto (S):", 1000),
            ("OD_recto", "Diámetro tramo recto (S):", od_sug),
            ("L_deriv", "Longitud tramo derivación (D):", 1000),
            ("OD_deriv", "Diámetro tramo derivación (D):", od_sug),
        ]
        if grado == 0:
            # la Tee es la primerísima pieza de la red: hace falta además
            # el tramo que la alimenta, porque no hay ningún tramo previo
            # del que "heredar" la entrada.
            campos_tee = [
                ("L_in", "Longitud tramo de entrada:", 1000),
                ("OD_in", "Diámetro tramo de entrada:", od_sug),
            ] + campos_tee

        datos = self._pedir_valores("Agregar Tee", campos_tee)
        if datos is None:
            return
        try:
            v = {k: float(x) for k, x in datos.items()}
        except ValueError:
            messagebox.showerror("Valor inválido", "Todos los campos deben ser números.")
            return

        usados = set(self._ids_actuales())
        if grado == 0:
            tee_id = self._siguiente_id_excluyendo(usados)
            usados.add(tee_id)
        else:
            tee_id = origen  # convierte el extremo libre existente en la Tee
        nodo_s = self._siguiente_id_excluyendo(usados)
        usados.add(nodo_s)
        nodo_d = self._siguiente_id_excluyendo(usados)

        decl_tee = (
            f'{tee_id}["Tee:OD2={_fmt(v["OD2"])},Alpha={_fmt(v["Alpha"])},'
            f'H1={_fmt(v["H1"])},H2={_fmt(v["H2"])},H3={_fmt(v["H3"])},Giro={_fmt(v["Giro"])}"]'
        )
        decl_s = f"{nodo_s}((Libre))"
        decl_d = f"{nodo_d}((Libre))"
        tramo_s = f'{tee_id} -->|L={_fmt(v["L_recto"])},OD={_fmt(v["OD_recto"])},S| {nodo_s}'
        tramo_d = f'{tee_id} -->|L={_fmt(v["L_deriv"])},OD={_fmt(v["OD_deriv"])},D| {nodo_d}'

        if grado == 0:
            tramo_in = f'{origen} -->|L={_fmt(v["L_in"])},OD={_fmt(v["OD_in"])}| {tee_id}'
            self._insertar_lineas([decl_tee, decl_s, decl_d], [tramo_in, tramo_s, tramo_d])
        else:
            self._insertar_lineas([decl_s, decl_d], [tramo_s, tramo_d], reemplazo_nodo=(tee_id, decl_tee))

        self._fijar_nodo_origen(nodo_s)
        self._dibujar_diagrama()
        self._escribir_log(
            f"Tee en {tee_id}: recto -> {nodo_s} / derivación -> {nodo_d}. "
            f"Nodo de origen actual: {nodo_s} (cambialo al combo si querés seguir por la derivación {nodo_d})."
        )

    def _agregar_codo(self):
        origen = self._origen_o_aviso()
        if origen is None:
            return
        grado = self._grado_nodo(origen)
        if grado >= 1 and not self._tramo_entrante(origen):
            self._advertir_nodo_no_disponible(origen)
            return

        od_sug = self._od_sugerido()
        campos = [
            ("R", "R (radio de curvatura):", 225),
            ("Angulo", "Ángulo del codo (°):", 90),
            ("Giro", "Giro (°):", 0),
            ("L_out", "Longitud tramo de salida:", 1000),
            ("OD_out", "Diámetro tramo de salida:", od_sug),
        ]
        if grado == 0:
            campos = [
                ("L_in", "Longitud tramo de entrada:", 1000),
                ("OD_in", "Diámetro tramo de entrada:", od_sug),
            ] + campos

        datos = self._pedir_valores("Agregar Codo", campos)
        if datos is None:
            return
        try:
            v = {k: float(x) for k, x in datos.items()}
        except ValueError:
            messagebox.showerror("Valor inválido", "Todos los campos deben ser números.")
            return

        usados = set(self._ids_actuales())
        if grado == 0:
            codo_id = self._siguiente_id_excluyendo(usados)
            usados.add(codo_id)
        else:
            codo_id = origen  # convierte el extremo libre existente en el Codo
        nodo_out = self._siguiente_id_excluyendo(usados)

        decl_codo = f'{codo_id}["Codo:R={_fmt(v["R"])},Angulo={_fmt(v["Angulo"])},Giro={_fmt(v["Giro"])}"]'
        decl_out = f"{nodo_out}((Libre))"
        tramo_out = f'{codo_id} -->|L={_fmt(v["L_out"])},OD={_fmt(v["OD_out"])}| {nodo_out}'

        if grado == 0:
            tramo_in = f'{origen} -->|L={_fmt(v["L_in"])},OD={_fmt(v["OD_in"])}| {codo_id}'
            self._insertar_lineas([decl_codo, decl_out], [tramo_in, tramo_out])
        else:
            self._insertar_lineas([decl_out], [tramo_out], reemplazo_nodo=(codo_id, decl_codo))

        self._fijar_nodo_origen(nodo_out)
        self._dibujar_diagrama()
        self._escribir_log(f"Codo en {codo_id} -> {nodo_out}.")

    def _agregar_reduccion(self):
        origen = self._origen_o_aviso()
        if origen is None:
            return
        grado = self._grado_nodo(origen)
        if grado >= 1 and not self._tramo_entrante(origen):
            self._advertir_nodo_no_disponible(origen)
            return

        od_sug = self._od_sugerido()
        # OD1/OD2 de la reducción fijan también el diámetro de los tramos
        # que entran y salen de ella (si no coincidieran no tendría sentido
        # geométrico: la reducción cambia justo de un diámetro al otro).
        campos = [
            ("OD1", "Diámetro de entrada (OD1):", od_sug),
            ("L", "Longitud de la reducción:", 250),
            ("OD2", "Diámetro de salida (OD2):", max(od_sug - 50, 50)),
            ("L_out", "Longitud tramo de salida:", 1000),
        ]
        if grado == 0:
            campos = [("L_in", "Longitud tramo de entrada:", 1000)] + campos

        datos = self._pedir_valores("Agregar Reducción", campos)
        if datos is None:
            return
        try:
            v = {k: float(x) for k, x in datos.items()}
        except ValueError:
            messagebox.showerror("Valor inválido", "Todos los campos deben ser números.")
            return

        usados = set(self._ids_actuales())
        if grado == 0:
            red_id = self._siguiente_id_excluyendo(usados)
            usados.add(red_id)
        else:
            red_id = origen  # convierte el extremo libre existente en la Reducción
        nodo_out = self._siguiente_id_excluyendo(usados)

        decl_red = f'{red_id}["Reduccion:OD1={_fmt(v["OD1"])},OD2={_fmt(v["OD2"])},L={_fmt(v["L"])}"]'
        decl_out = f"{nodo_out}((Libre))"
        tramo_out = f'{red_id} -->|L={_fmt(v["L_out"])},OD={_fmt(v["OD2"])}| {nodo_out}'

        if grado == 0:
            tramo_in = f'{origen} -->|L={_fmt(v["L_in"])},OD={_fmt(v["OD1"])}| {red_id}'
            self._insertar_lineas([decl_red, decl_out], [tramo_in, tramo_out])
        else:
            self._insertar_lineas([decl_out], [tramo_out], reemplazo_nodo=(red_id, decl_red))

        self._fijar_nodo_origen(nodo_out)
        self._dibujar_diagrama()
        self._escribir_log(f"Reducción en {red_id} -> {nodo_out}.")


def main():
    archivo_inicial = sys.argv[1] if len(sys.argv) > 1 else None
    app = App(archivo_inicial)
    app.mainloop()


if __name__ == "__main__":
    main()
