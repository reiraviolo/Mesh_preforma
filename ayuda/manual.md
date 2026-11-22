# Manual: malla de derivación (tee) — media pieza + espejado

Este manual explica, paso a paso, qué hace el notebook `malla_derivacion.ipynb` y por qué.
Está pensado para ir aprendiendo: leé un paso, mirá la celda correspondiente, cambiá un parámetro y volvé a ejecutar.

---

## 0. Qué problema resuelve

Queremos una **malla de cuadriláteros** (elementos de cáscara/shell) para una **tee**: un caño principal con una derivación que sale perpendicular.

Para no mallar todo, se aprovecha la **simetría**: se malla solo la mitad (y ≥ 0) y después se **espeja**.

### Geometría

| Elemento | Eje | Radio |
|---|---|---|
| Caño principal (*run*) | X | `R` |
| Derivación (*branch*) | Z | `r` (con `r < R`) |

- Plano de simetría: **zx** (y = 0).
- La derivación arranca en la **curva de intersección** con el caño principal y sube hasta `z = H`.

> Nota: en este notebook solo se malla la **derivación**. El caño principal se dibuja como superficie de referencia, pero todavía no tiene malla.

---

## 1. Imports y parámetros (celda 1)

```python
R = 100.0; r = 50.0; Lx = 250.0; H = 250.0
n_theta = 16; n_z = 10; razon_z = 1.20
```

**Geometría**
- `R`: radio del caño principal.
- `r`: radio de la derivación (tiene que ser menor que `R`).
- `Lx`: semilongitud del caño principal (va de −Lx a +Lx). Solo afecta al dibujo.
- `H`: cota Z del extremo libre de la derivación.

**Malla**
- `n_theta`: cantidad de elementos en la media circunferencia (θ de 0 a π).
- `n_z`: cantidad de elementos a lo largo de la derivación.
- `razon_z`: cuánto más grande es cada elemento respecto del anterior (1.20 = crece 20 % por elemento; 1 = todos iguales).

Los `assert` al final frenan el programa si la geometría no tiene sentido (`r ≥ R` o `H` por debajo del caño).

**Para probar:** poné `razon_z = 1.0` y fijate que los elementos pasan a ser todos iguales.

---

## 2. Curva de intersección (celda 2)

Es la línea donde se cruzan los dos cilindros. Se obtiene así:

- La derivación es un cilindro de radio `r` con eje Z, entonces un punto de su superficie es
  `x = r·cosθ`, `y = r·sinθ`.
- El caño principal cumple `y² + z² = R²`, entonces `z = √(R² − y²)`.

Reemplazando `y = r·sinθ` queda `z = √(R² − (r·sinθ)²)`.

La función `curva_interseccion(R, r, theta)` devuelve un arreglo `(n, 3)` con `(x, y, z)`.

Con los valores por defecto:
- en θ = 0 (y = 0): `z = 100` (la cresta del caño principal);
- en θ = 90° (y = 50): `z ≈ 86.6` (la intersección "baja" a los costados).

Esa diferencia de altura es la que hace que la base de la derivación no sea plana.

---

## 3. Distribuciones de nodos (celda 3)

Dos funciones chicas que deciden **dónde** van los nodos.

**`distribucion_theta(n)`**
Devuelve `n+1` ángulos uniformes entre 0 y π. Si más adelante querés refinar la malla en alguna zona (por ejemplo cerca de θ = 0, donde está la cresta), este es el lugar para modificar.

**`distribucion_geometrica(n, razon)`**
Devuelve `n+1` valores de 0 a 1 donde cada intervalo es `razon` veces el anterior. La fórmula es la suma de una serie geométrica:

```
s_k = (razon^k − 1) / (razon^n − 1)
```

Con `n = 10` y `razon = 1.2` los valores son aproximadamente:

```
0, 0.039, 0.085, 0.140, 0.207, 0.287, 0.383, 0.498, 0.636, 0.801, 1.0
```

El primer elemento mide ~4 % del largo y el último ~20 %: **elementos chicos cerca de la intersección** (donde hay concentración de tensiones) y grandes hacia el extremo libre.

---

## 4. Malla de la derivación (celda 4)

Es el corazón del programa: `malla_derivacion_mitad(...)`.

### 4.1 Grilla estructurada

Se arma una grilla de `(n_theta+1) × (n_z+1)` nodos:
- índice `i` → ángulo θ (0 … n_theta)
- índice `j` → altura (0 … n_z), con `j = 0` en la intersección y `j = n_z` en el extremo libre

Numeración de nodos:

```
id(i, j) = i·(n_z + 1) + j
```

### 4.2 Posición de cada nodo

Para cada ángulo θᵢ se calcula primero la altura de la intersección `z_int[i]`. Después el nodo `j` queda en:

```
x = r·cosθ
y = r·sinθ
z = z_int + (H − z_int)·s_j
```

O sea: cada "meridiano" va desde la curva de intersección hasta `z = H`, repartiendo los nodos con la distribución geométrica `s`.

La línea que dice `np.where(abs(y) < 1e-12·r, 0, y)` limpia el error numérico: `sin(π)` no da exactamente 0, y necesitamos que los nodos del plano de simetría tengan **y = 0 exacto** para poder detectarlos después al espejar.

### 4.3 Elementos (cuadriláteros)

Cada celda de la grilla es un cuadrilátero con 4 nodos:

```
(i, j) → (i+1, j) → (i+1, j+1) → (i, j+1)
```

Ese orden es **antihorario visto desde afuera**, es decir, la normal del elemento apunta hacia afuera del caño. (Se puede verificar: la dirección de θ creciente × la dirección Z da la dirección radial saliente.) Mantener la normal consistente es importante para cargas de presión en un modelo de elementos finitos.

### 4.4 Grupos de nodos

La función también devuelve un diccionario `grupos` con índices de nodos útiles:

- `interseccion`: nodos de la base (j = 0). Ahí se va a pegar la malla del caño principal.
- `extremo`: nodos del borde libre (j = n_z). Ahí se aplicarían cargas o condiciones de borde.
- `simetria`: nodos con y = 0 (los de i = 0 e i = n_theta). Ahí va la condición de simetría.

**Con los valores por defecto:** 17 × 11 = **187 nodos** y 16 × 10 = **160 elementos**.

---

## 5. Espejado (celda 5)

`espejar_malla(nodos, quads, valores, cols_y)` completa la pieza reflejando respecto del plano y = 0.

Paso a paso:

1. **Detectar nodos en el plano** (`|y| < tol`). Esos **no se duplican**, porque el espejo de un nodo con y = 0 es él mismo. Si se duplicaran quedaría una costura abierta.
2. **Crear los nodos nuevos**: se copian los demás cambiando `y → −y`.
3. **Mapa de espejo** (`mapa`): para cada nodo original `k`, `mapa[k]` dice cuál es su espejo en la malla completa. Para los nodos del plano, `mapa[k] = k`.
4. **Espejar los elementos**: se reemplazan los índices usando `mapa` y se **invierte el orden** (`[:, ::-1]`). Esto es necesario porque un espejo invierte la orientación; sin esto las normales de la mitad nueva apuntarían hacia adentro.
5. **Espejar valores nodales** (opcional): los valores se copian tal cual, salvo las columnas que sean la **componente Y de un vector** (desplazamiento, normal, fuerza…), que cambian de signo. Esas columnas se indican con `cols_y`.

**Con los valores por defecto:** 22 nodos en el plano, se agregan 165 → **352 nodos** y **320 elementos**.

---

## 6. Gráfico (celda 6)

`graficar(...)` dibuja con matplotlib 3D:

- el caño principal (azul tenue);
- la derivación completa (naranja tenue, solo de referencia);
- la curva de intersección completa en rojo;
- la malla de la derivación (cuadriláteros con borde negro);
- opcionalmente, la mitad espejada en celeste (`mostrar_espejo=True`).

Si querés rotar con el mouse, usá `%matplotlib widget` (requiere `ipympl`).

---

## 7. Ejecución (celda 7)

```python
nodos, quads, grupos = malla_derivacion_mitad(R, r, H, n_theta, n_z, razon_z)
graficar(R, r, Lx, H, nodos, quads)
```

Genera la media malla, imprime cantidad de nodos/elementos y la dibuja.

---

## 8. Espejar y ver la pieza completa (celda 8, opcional)

Crea un campo nodal de prueba con 4 columnas: `z`, `x`, `y`, `z`. La columna de índice 2 es la componente Y de un vector, por eso se pasa `cols_y=(2,)`. Después llama a `espejar_malla` y grafica la pieza completa.

---

## Orden de ejecución

La primera vez, de arriba hacia abajo: 1 → 2 → 3 → 4 → 5 → 6 → 7 → 8.
Después se puede cambiar cualquier celda de función y volver a ejecutar solo esa y las de ejecución (7 y 8).

## Ejercicios sugeridos

1. Cambiá `n_theta` a 24 y `n_z` a 15. ¿Cuántos nodos y elementos da? (Fórmula: nodos = (n_theta+1)(n_z+1), elementos = n_theta·n_z.)
2. Poné `razon_z = 1.0` y compará con `1.3`.
3. Poné `r = 90` (casi igual a `R`) y mirá cómo se deforma la curva de intersección.
4. En `distribucion_theta`, concentrá nodos cerca de θ = 0 (por ejemplo con `np.pi * np.linspace(0,1,n+1)**1.5`).

## Qué falta (próximos pasos del proyecto)

- Mallar el **caño principal** (con el agujero de la intersección) de forma que sus nodos coincidan con los de `grupos['interseccion']`.
- Después, repetir la idea para **codos, tubos y reducciones**, que son las otras piezas del proyecto.

---

> El notebook fue el prototipo de la tee. El repo `Mesh_preforma` ya implementa el programa completo (tubo, codo, reducción, tee, ensamblado y exportación). La **Parte 2** explica ese programa real, empezando por lo más simple.

---

# Parte 2 — El programa `Mesh_preforma`

## 9. Visión general: qué hace el programa de punta a punta

El flujo está en `main.py` y se ejecuta con:

```
uv run python main.py ejemplos/Caso_1.md 48
```

Pasos (cada uno es una función de un módulo de `src/`):

| # | Paso | Módulo | Qué hace |
|---|---|---|---|
| 1 | `leer_mermaid` | `mermaid_reader.py` | Lee el diagrama Mermaid: qué accesorios hay y cómo se conectan. |
| 2 | `ensamblar` | `ensamblador.py` | Genera la malla de cada pieza y la ubica en el espacio. |
| 3 | `calcular_posiciones` | `red.py` | Calcula dónde está cada nodo de la red y hacia dónde apunta cada tramo. |
| 4 | `fusionar_piezas` | `ensamblador.py` | Une todas las mallas en una sola, unificando nodos coincidentes. Devuelve también las **zonas de soldadura**. |
| 5 | `escribir_malla` | `exportador.py` | Escribe el archivo `.malla`. |
| 6 | `graficar_*` | `graficador.py` | Muestra el resultado en 3D. |

La idea central: **cada pieza se mallea una sola vez en su propio sistema de coordenadas**, y después se rota y traslada a su lugar. Es mucho más simple que mallar toda la red de golpe.

### El archivo de entrada (Mermaid)

Ejemplo `Caso_1.md`:

```
1((Libre))
2["Tee:OD2=100,Alpha=60,H1=200,H2=200,H3=211,Giro=30"]
3["Codo:R=225,Angulo=45,Giro=45"]
...
1 -->|L=1000,OD=150| 2
2 -->|L=2000,OD=150,S| 3
```

- Los **nodos** (`2["Tee:..."]`) son accesorios o extremos libres.
- Las **flechas** (`1 --> 2`) son tramos de caño recto, con su largo `L` y diámetro `OD`.
- `S` / `D` en un tramo que sale de una Tee indica si es la salida directa (*Straight*) o la derivación (*Derivación*).

Es decir: los accesorios son los nodos del grafo y los **tubos son las aristas**.

## 10. El tubo (`src/parts/tubo.py`)

Es la pieza más simple y la base de las demás. La función pública es `generar_malla_tubo(OD, L, N)`; el trabajo lo hace `_malla_cilindro`.

**Entradas:** `OD` (diámetro exterior), `L` (largo), `N` (divisiones del perímetro).

### 10.1 Tamaño de los elementos

```python
delta = pi * D / N                       # lado del elemento sobre el perímetro
N_segmentos_L = max(1, round(L / delta)) # cuántos elementos caben a lo largo
paso_L = L / N_segmentos_L               # paso axial real
```

Se busca que los elementos sean **casi cuadrados**: el lado circunferencial es `π·D/N`, y a lo largo se elige la cantidad de segmentos que más se acerque a ese tamaño. Después se ajusta `paso_L` para que la suma dé exactamente `L`.

Ejemplo (`OD=200, L=1000, N=40`): `delta = 15.7`, `L/delta = 63.7` → **64 segmentos**, `paso_L = 15.625`.

### 10.2 Nodos

Los nodos forman una matriz de `N_filas × N` (anillos × puntos por anillo), numerados por fila:

```
id = fila * N + columna
```

Para cada nodo:

```
angulo = columna * 2π / N
x = x0 + fila * paso_L        (el eje del caño es X)
y = (D/2) * cos(angulo)
z = (D/2) * sin(angulo)
```

Hay `N_filas = N_segmentos_L + 1` anillos.

### 10.3 Elementos y cierre del anillo

Cada cuadrilátero toma 4 nodos de una "ventana" 2×2 de la matriz. El truco está en la matriz auxiliar `buffer`, que tiene **una columna extra** (`N+1`) que repite la columna 0:

```python
buffer[:, :N] = matriz
buffer[:, N]  = matriz[:, 0]   # cierra el anillo
```

Así el último elemento de cada anillo se conecta con el primero sin duplicar nodos.

Al final, `[:, [0, 1, 3, 2]]` reordena los 4 nodos para que den la vuelta al elemento en orden (perimetral) en lugar del orden de lectura de la ventana.

### 10.4 Orientación

El orden resultante tiene la normal **hacia afuera** del caño (lo verifiqué con el primer elemento: normal ≈ dirección radial). Es la misma convención que usaba el notebook.

### 10.5 Resultado

Devuelve un diccionario con `nodos`, `elementos`, `grupo` (ceros, un solo grupo "tubo") y `divisiones` (`OD, L, N, N_filas, paso_L`).

Con `OD=200, L=1000, N=40`: **2600 nodos** y **2560 elementos** (40 × 64).

> Por qué guarda `N_filas` y `paso_L`: los accesorios vecinos necesitan saber cuántos nodos tiene el anillo de unión (`N`) para que las mallas **coincidan** al fusionarlas.

## 11. Poner cada pieza en su lugar (`src/transformaciones.py`)

Todas las piezas nacen con la entrada apuntando a un eje local conocido: **+X** para las cuatro (tubo, codo, reducción y tee). Hay que rotarlas para que apunten donde manda la red.

### 11.1 `rotar_eje` — rotación de Rodrigues

Rota un vector `v` un ángulo θ alrededor de un eje cualquiera:

```
v_rot = v·cosθ + (eje × v)·sinθ + eje·(eje·v)·(1 − cosθ)
```

Idea: se separa `v` en la parte **paralela** al eje (no cambia) y la **perpendicular** (gira). La usa `red.py` para calcular hacia dónde apunta cada tramo según `Alpha`, `Angulo` y `Giro`.

### 11.2 `_matriz_rotacion` — la misma fórmula, en forma de matriz 3×3

Construye `R = I + sinθ·K + (1−cosθ)·K²`, donde `K` es la matriz del producto cruz con el eje. Sirve para rotar **todos los nodos de una vez** con una multiplicación de matrices.

### 11.3 `alinear_pieza` — orientar una pieza completa

Recibe dos pares de ejes:
- **local:** hacia dónde apunta la entrada de la pieza y su eje normal (el perpendicular al plano donde se desvía; en el ensamblador se usa **−Y** para el codo y la tee, ver sección 28).
- **global:** hacia dónde debe apuntar según `red.py`.

Hace la rotación en **dos pasos**:
1. **Rotación mínima** que lleva la entrada local a la entrada global (eje de giro = producto cruz entre ambas).
2. **Giro (twist)** alrededor de la entrada ya alineada, para que el eje normal también coincida. Es lo que implementa el parámetro `Giro`.

Casos especiales que contempla: si ya están alineadas no rota; si son exactamente opuestas (180°) elige cualquier eje perpendicular.

Para tubo y reducción alcanza el paso 1, porque la sección es circular y el giro no cambia nada.

### 11.4 `traslacion` y `escalado`

Helpers triviales: sumar un vector a todos los nodos o multiplicarlos por un factor.

---

---

# Parte 3 — Codo, reducción y lectura del Mermaid

## 12. El codo (`src/parts/codo.py`)

Un codo es un **tubo que se curva**: en vez de avanzar en línea recta, el centro de cada anillo recorre un arco de circunferencia (un toro). La función pública es `generar_malla_codo(OD, R, angulo, N)`.

**Entradas:** `OD`, `R` (radio de curvatura, medido al **eje** del caño), `angulo` (grados, por defecto 90) y `N`.

### 12.1 Cuántos anillos

```python
theta = deg2rad(angulo)
L_arco = R * theta                               # largo del arco sobre el eje
N_filas_L = max(1, int(L_arco // (pi*D/N)))      # cuántos elementos a lo largo
N_filas = N_filas_L + 1
```

Es la misma idea que en el tubo (elementos casi cuadrados), pero con dos diferencias:
- el largo es el **arco** `R·θ`, no `L`;
- usa división entera `//` (redondea hacia abajo), mientras que el tubo usa `round`.

Ejemplo (`OD=200, R=300, angulo=90, N=25`): `L_arco = 471.2`, `π·D/N = 25.1` → **18 elementos**, 19 anillos.

### 12.2 Elementos y cierre del anillo

Idéntico al tubo: matriz `N_filas × N`, columna extra en `buffer` para cerrar el anillo, y un reordenamiento final. Acá el reordenamiento es `[2, 3, 1, 0]` y no `[0, 1, 3, 2]`. El comentario del código explica por qué: con la parametrización actual ese orden recorre el perímetro al revés que el tubo, y se corrigió porque `control_malla` marcaba **48 elementos con orientación inconsistente** en la red ensamblada; al invertirlo pasó a 0. (Lo cual refuerza la idea de que **todas las piezas tienen que tener la misma orientación**, o la malla fusionada queda con normales mezcladas.)

### 12.3 Nodos: cómo se construye cada punto

Para cada nodo hay dos parámetros:
- `phi`: ángulo alrededor de la sección (0 … 2π), como en el tubo;
- `t`: ángulo recorrido a lo largo del arco (0 … θ).

Con ellos se calcula:

```python
centro = [R*sin(t), 0, R*(1 - cos(t))]      # punto del eje en ese anillo
e1     = [-sin(t), 0, cos(t)]               # primer vector de la sección
e2     = [0, 1, 0]                          # segundo vector de la sección
nodo   = centro + (D/2) * (cos(phi)*e1 + sin(phi)*e2)
```

Interpretación:
- **`centro`**: recorre un arco de radio `R` en el plano X-Z. Arranca en el origen y su tangente inicial es **+X** (igual que el tubo y la reducción, por eso todos acoplan sin rotaciones raras). Para 90° termina en `(R, 0, R)`: lo verifiqué, con `R=300` el último anillo queda centrado en `(300, 0, 300)`.
- **`e1` y `e2`**: dos vectores perpendiculares a la tangente que arman la circunferencia de la sección. `e1` **gira junto con la curva** (queda siempre en el plano X-Z y perpendicular a la tangente); `e2 = +Y` es fijo, porque es perpendicular al plano de giro.
- La sección es el círculo `centro + (D/2)(cosφ·e1 + sinφ·e2)`.

### 12.4 Intradós y extradós: los elementos no son iguales

En un codo, los elementos del lado **interno** (intradós, el de adentro de la curva) quedan más angostos en la dirección del arco que los del lado **externo** (extradós). Con `R=300` y `OD=200`, el avance entre anillos mide **17.45** en el intradós contra **34.83** en el extradós, una relación de ~2, que coincide con `(R + D/2) / (R − D/2) = 400/200`.

Los elementos casi cuadrados del tubo acá son un compromiso: el tamaño se calcula con el eje (`R·θ`), así que en el eje son casi cuadrados, en el intradós quedan más finos y en el extradós más alargados. Es normal en mallas de codo; conviene tenerlo presente si se usan relaciones de radio chicas.

### 12.5 Resultado

Resultado con los valores de arriba: **475 nodos** (25 × 19) y **450 elementos** (25 × 18).

## 13. La reducción (`src/parts/reduccion.py`)

Es un tubo cuyo diámetro cambia de `OD_1` (entrada, x = 0) a `OD_2` (salida, x = L). Función pública: `generar_malla_reduccion(OD_1, OD_2, L, N)`.

### 13.1 Tres tramos

La reducción **no** es un cono puro: tiene tres tramos.

```
recto con OD_1  →  cono (OD_1 a OD_2)  →  recto con OD_2
   L_D1 = 0.4·OD_1       L_cono               L_D2 = 0.4·OD_2
```

- Los tramos rectos miden 0.4 veces su diámetro. Sirven para que la unión con el caño de al lado quede sobre una zona cilíndrica, sin la curvatura del cono.
- Si los dos tramos rectos suman más de la mitad de `L`, se **achican proporcionalmente** (`factor = 0.5·L / (L_D1 + L_D2)`) para que el cono siempre tenga al menos la mitad del largo.
- El cono ocupa lo que queda: `L_cono = L − L_D1 − L_D2`.

### 13.2 Divisiones por tramo

Para cada tramo se usa el mismo criterio de elementos casi cuadrados, con el diámetro propio del tramo (el del cono es el promedio de los dos):

```python
def divisiones_tramo(largo, D):
    delta = pi * D / N
    return max(1, int(largo // delta))
```

Ejemplo del `Caso_1` (`OD1=150, OD2=100, L=250, N=48`): `L_D1 = 60`, `L_cono = 150`, `L_D2 = 40`, con **6 + 18 + 6** divisiones.

### 13.3 Posición de los anillos y diámetro de cada uno

Se arma el vector de posiciones `x_anillos` concatenando los tres tramos (con `[1:]` para no repetir el anillo compartido) y el diámetro de cada anillo se obtiene **interpolando**:

```python
D_anillos = np.interp(x_anillos, [0, x1, x2, L], [OD_1, OD_1, OD_2, OD_2])
```

Eso da diámetro constante en los rectos y variación lineal en el cono.

### 13.4 Nodos y elementos

Igual que el tubo, pero con el diámetro de cada fila (`D_actual`) en vez de uno solo: `y = (D_actual/2)·cos(ángulo)`, `z = (D_actual/2)·sin(ángulo)`. Los elementos usan el mismo orden `[0, 1, 3, 2]` y el mismo cierre de anillo.

Ejemplo (`OD_1=100, OD_2=50, L=150, N=25`): 16 anillos, **400 nodos**, **375 elementos**.

> Dato práctico: como la entrada es siempre **+X**, tubo, reducción y tee (y ahora el codo) empiezan por el mismo lado, y por eso `alinear_pieza` solo tiene que rotarlas a su orientación final.

## 14. El lector de Mermaid (`src/mermaid_reader.py`)

Convierte el texto del diagrama en dos estructuras de datos. **No genera geometría ni mallas**; solo interpreta texto.

### 14.1 Tres tipos de línea, tres expresiones regulares

| Tipo de línea | Ejemplo | Patrón |
|---|---|---|
| Nodo libre | `1((Libre))` | `PATRON_NODO_LIBRE` |
| Nodo con parámetros | `2["Tee:OD2=100,Alpha=60,..."]` | `PATRON_NODO_PARAMS` |
| Tramo | `1 -->|L=1000,OD=150| 2` | `PATRON_TRAMO` |

Se ignoran las líneas que empiezan con `` ``` ``, `flowchart` o `graph`. Si una línea no coincide con ningún patrón, **se avisa por pantalla** y se ignora (no falla en silencio).

### 14.2 Parámetros: `_parsear_params`

Convierte `"OD2=100,Alpha=60,S"` en un diccionario:

```python
{"OD2": 100, "Alpha": 60, "S": True}
```

Reglas:
- `clave=valor` → número si se puede (entero si no tiene decimales, si no `float`), o texto si no es numérico.
- También acepta `clave:valor` (por si alguien usa `:` en vez de `=`).
- Un token **sin** `=` ni `:` (como `S` o `D`) se guarda como bandera `True`.

### 14.3 Resultado: `(nodos, tramos)`

Para el `Caso_1.md`:

```python
nodos  = {"1": {"tipo": "Libre", "params": {}},
          "2": {"tipo": "Tee", "params": {"OD2": 100, "Alpha": 60, "H1": 200, "H2": 200, "H3": 211, "Giro": 30}},
          "3": {"tipo": "Codo", "params": {"R": 225, "Angulo": 45, "Giro": 45}}, ...}

tramos = [{"origen": "1", "destino": "2", "params": {"L": 1000, "OD": 150}},
          {"origen": "2", "destino": "3", "params": {"L": 2000, "OD": 150, "S": True}}, ...]
```

Un nodo que aparece en un tramo pero nunca se definió se agrega como `Libre` (con aviso).

### 14.4 Dos funciones

- `parsear_mermaid(texto)`: trabaja sobre texto en memoria. La usa la interfaz gráfica, que tiene el texto en un editor y no un archivo.
- `leer_mermaid(path)`: abre el archivo y llama a la anterior. La usa `main.py`.

---

---

# Parte 4 — La geometría de la red (`src/red.py`)

Hasta ahora tenemos piezas sueltas (Parte 2 y 3) y un grafo en texto (sección 14). `red.py` es el puente: **recorre el grafo y calcula, para cada nodo, dónde está y hacia dónde apunta el caño que llega a él**. Todavía no genera mallas; es un paso intermedio para poder verificar el ruteo (con prints o graficando los ejes) antes de mallar.

La función pública es `calcular_posiciones(nodos, tramos)`.

## 15. Qué devuelve

Un diccionario `{id: {"pos", "dir_in", "ref_in"}}`:

| Campo | Significado |
|---|---|
| `pos` | Posición 3D del nodo (centro del accesorio). |
| `dir_in` | Dirección del tramo que **entra** al nodo (vector unitario). |
| `ref_in` | Vector perpendicular a `dir_in` que define **dónde está Giro = 0** para el accesorio de ese nodo (ver sección 18). |

La raíz tiene `dir_in = None` y `ref_in = None`.

Resultado real para `Caso_1.md` (redondeado):

| Nodo | Tipo | `pos` | `dir_in` |
|---|---|---|---|
| 1 | Libre | (0, 0, 0) | — |
| 2 | Tee | (1000, 0, 0) | (1, 0, 0) |
| 3 | Codo | (3000, 0, 0) | (1, 0, 0) |
| 6 | Libre (derivación) | (1500, 750, 433) | (0.5, 0.75, 0.433) |
| 4 | Reducción | (3707, 183, 683) | (0.707, 0.183, 0.683) |
| 5 | Libre | (4414, 366, 1366) | (0.707, 0.183, 0.683) |

## 16. El algoritmo, paso a paso

### Paso 1 — Validar la topología (`_validar_grados`)

Cuenta cuántos tramos toca cada nodo y lo compara con lo que corresponde a su tipo:

```python
GRADOS_ESPERADOS = {"Tee": 3, "Codo": 2, "Reduccion": 2, "Libre": 1}
```

Una tee tiene 3 bocas, un codo y una reducción 2 y un extremo libre 1. Si algo no cierra, se lanza un error claro; por ejemplo, un codo con un solo caño da `Nodo 2 (Codo): tiene 1 tramo(s), se esperaba(n) 2.`

### Paso 2 — Buscar la raíz

Se arma, para cada nodo, la lista de tramos que **salen** de él, y el conjunto de nodos que tienen tramo **entrante**. La raíz es el único nodo que no tiene entrante. Si hay cero o más de una, error: la red tiene que ser un **árbol con una sola raíz** (sin lazos cerrados).

### Paso 3 — Recorrer el árbol

Es un recorrido en profundidad con una pila (`pendientes`):

1. La raíz se coloca en `(0, 0, 0)`.
2. Se saca un nodo de la pila y, por cada tramo que sale de él:
   - se calcula la **dirección de salida** (y el ref que se propaga);
   - la posición del nodo destino es `pos_actual + dir_salida · L`;
   - se guarda el destino y se lo agrega a la pila.
3. Se repite hasta vaciar la pila.

Detalle importante: **`L` es la distancia entre los centros de los nodos** a lo largo del eje. El primer tramo sale de la raíz con dirección fija `(1, 0, 0)` por convención (parámetro `direccion_inicial`). Cómo usa `ensamblador.py` ese `L` para generar los tubos (acortándolo con los "recortes" de los accesorios vecinos) se explica en la sección 29.

## 17. Cómo cambia la dirección en cada accesorio (`_direccion_salida`)

Cada tipo de nodo tiene su regla:

| Tipo | Dirección de salida | Ref que se propaga |
|---|---|---|
| **Reducción** | La misma que entra (solo cambia el diámetro). | `ref_in` sin cambios. |
| **Codo** | `dir_in` girada `Angulo` grados alrededor de un **eje de desvío**. | El propio eje de desvío. |
| **Tee, salida `S`** | La misma que entra (el tronco sigue derecho). | El eje de desvío (ver abajo). |
| **Tee, salida `D`** | `dir_in` girada `Alpha` grados alrededor del eje de desvío. | El eje de desvío. |
| **Libre** | No tiene salida: es una hoja o la raíz. | — |

Una tee con un tramo saliente sin flag `S` ni `D` produce un error: no se puede adivinar cuál es el tronco.

### Por qué la tee propaga el eje en las dos salidas

Aunque el tronco `S` sigue derecho, la tee es **un solo cuerpo rígido** (tronco + derivación se generan y rotan juntos). El `Giro` de la tee rota todo el cuerpo, incluida la cara de salida del tronco. Si el tramo recto siguiente usara el `ref_in` sin rotar, sus nodos no coincidirían nodo a nodo con la cara real del tronco. Por eso se propaga el eje ya rotado.

## 18. El parámetro `Giro`: eje de desvío y referencia

Esta es la parte más sutil, así que va despacio.

### El problema

Un codo de `Angulo = 90` cambia la dirección en 90°, pero **¿hacia qué lado**? Hay infinitas direcciones a 90° de `dir_in` (todas las del plano perpendicular). `Giro` elige una.

### La solución: un eje de desvío

La salida se obtiene rotando `dir_in` alrededor de un **eje de desvío** perpendicular a `dir_in`. Ese eje se calcula con `_eje_desvio`:

```python
eje = rotar_eje(ref_in, dir_in, giro_grados)
```

Es decir: se toma `ref_in` (el eje que corresponde a `Giro = 0`) y se lo **gira `Giro` grados alrededor de `dir_in`** con la fórmula de Rodrigues (sección 11.1). Después, `dir_in` se gira `Angulo` (o `Alpha`) alrededor de ese eje.

Este eje queda perpendicular tanto a la dirección que entra como a la que sale. Es la **normal al plano de desvío**.

### `Giro` no es respecto de una horizontal global

`Giro` se mide respecto del **plano de la pieza anterior** (su propio plano de desvío, o el que se venía arrastrando si la pieza anterior era recta). Por eso cada nodo lleva `ref_in` heredado. Así, un `Giro = 0` en un codo significa "seguí doblando en el mismo plano que el codo anterior", y no "doblá en un plano horizontal fijo".

### El primer accesorio: `_ref_inicial`

El primer accesorio de la red no tiene pieza anterior de la que heredar. Para él se arma un `ref` por convención: el **eje Z global proyectado perpendicularmente a `dir_in`**. Si `dir_in` ya es casi vertical (casi paralela a Z), Z no sirve y se usa X global.

Con la dirección inicial `+X`, el `ref_inicial` es `(0, 0, 1)` (eje Z vertical). Como el eje de desvío es vertical, el codo dobla **en el plano horizontal**.

### Verificación numérica

Probé un codo de 90° entrando en `+X` con el ref inicial, variando `Giro`:

| `Giro` | Eje de desvío | Salida |
|---|---|---|
| 0° | (0, 0, 1) | **+Y** (horizontal) |
| 90° | (0, −1, 0) | **+Z** (vertical hacia arriba) |
| 180° | (0, 0, −1) | **−Y** (horizontal, lado opuesto) |

Es decir: `Giro = 0` da una salida horizontal y `Giro = 90` una vertical, con la regla de la mano derecha alrededor de `dir_in`.

### Ejemplo: la tee del `Caso_1`

`Tee: Alpha=60, Giro=30`, entrando en `+X` con `ref_in = (0, 0, 1)`:

1. Eje de desvío = `(0, 0, 1)` girado 30° alrededor de `+X` = `(0, −0.5, 0.866)`.
2. Dirección de la derivación (`D`) = `(1, 0, 0)` girada 60° alrededor de ese eje = `(0.5, 0.75, 0.433)`.
3. Posición del nodo 6 = `(1000, 0, 0) + 1000 · (0.5, 0.75, 0.433)` = `(1500, 750, 433)`. ✔ Coincide con la tabla de la sección 15.
4. El tronco `S` sigue en `+X`, con ref = el mismo eje `(0, −0.5, 0.866)`. Por eso el codo siguiente (nodo 3, `Giro = 45`) ya hereda esa referencia y no la original.

## 19. Cómo se conecta con lo que viene

- `red.py` entrega `dir_in` y los ejes de desvío a `ensamblador.py`, que los usa con `alinear_pieza` (sección 11.3) para orientar cada pieza (confirmado en la Parte 6, secciones 28 y 29).
- `main.py` también usa `posiciones` para que `graficador.graficar_ejes` dibuje la línea de eje de cada tramo sobre la malla.

---

---

# Parte 5 — La tee (`src/parts/tee.py`)

Es el archivo más grande (1224 líneas) porque es la pieza difícil: hay que mallar **dos cilindros que se cruzan**, con una derivación inclinada un ángulo cualquiera. La idea de fondo es **la misma del notebook** (Parte 1: media pieza + espejado), pero generalizada y con muchos pasos para mejorar la calidad de la malla.

La función pública es `generar_malla_tee(OD_1, OD_2, H_1, H_2, H_3, Alpha_2, ...)`. En `main.py`/`red.py` estos parámetros salen de la línea `Tee:OD2=..,Alpha=..,H1=..,H2=..,H3=..`.

## 20. Geometría y parámetros

| Parámetro | Significado |
|---|---|
| `OD_1` | Diámetro del caño **directo** (el tronco; eje X). |
| `OD_2` | Diámetro de la **derivación** (`OD_2 ≤ OD_1`). |
| `H_1` | Distancia de la intersección a la boca de **ingreso** del directo. |
| `H_3` | Distancia de la intersección a la boca de **salida** del directo. |
| `H_2` | Distancia de la intersección a la boca de la derivación, **sobre su eje**. |
| `Alpha_2` | Ángulo de la derivación respecto del eje X, en grados (entre 0 y 180). |

Diferencia con el notebook: allá la derivación era perpendicular (eje Z). Acá su eje es

```python
eje_s  = (cos α, 0, sin α)      # eje de la derivación, en el plano X-Z
perp_s = (-sin α, 0, cos α)     # perpendicular a eje_s en ese plano
```

Con `α = 90°` queda el caso del notebook.

**Parámetros de discretización** (los que más se usan):

- `N_directo` (80 por defecto): divisiones de la circunferencia del directo; se redondea a múltiplo de 4.
- `N_s`: divisiones de la derivación. Si es `None`, se elige para que los elementos sean ~cuadrados: `N_s = 2·round(π·R2 / tam)`.
- `razon`: crecimiento geométrico de los elementos a lo largo de la derivación (igual que en el notebook).
- `capas_base` (3): cantidad de capas de elementos en la **corona** (ver sección 22).
- `rel_capa` (0.75): relación de aspecto de esas capas.
- `ANGULO_VERTICE` (40°): giro de la curva a partir del cual se considera un **vértice**.
- `umbral_cuna`, `alinear`: pasos de mejora de calidad (pasos 6 y 7 de la sección 22).

## 21. Los 3 grupos de elementos

Cada elemento queda etiquetado en un arreglo `grupo`:

| Grupo | Qué es |
|---|---|
| **0** | Derivación. |
| **1** | **Corona**: las capas de elementos que rodean la curva de intersección **sobre el caño directo**. |
| **2** | Resto del caño directo. |

Esto permite después tratar cada zona distinto (por ejemplo, suavizar solo el grupo 2).

## 22. El proceso, paso a paso

Los pasos 1 y 1b son la preparación (la curva de intersección) y los pasos 2 a 8 son las siete etapas que ejecuta `generar_malla_tee` en el código. Como en el notebook, se trabaja con la **mitad y ≥ 0** y al final se espeja.

### Paso 1 — Curva de intersección (`punto_curva`)

Es la fórmula del notebook, generalizada. Si `φ` es el ángulo alrededor de la derivación:

```python
y = R2 · cos φ
z = √(R1² − y²)                      # el punto está sobre el caño directo
x = (z · cos α − R2 · sin φ) / sin α
```

Se deduce igualando un punto de la derivación `t·eje_s + R2·cosφ·ŷ + R2·sinφ·perp_s` con la ecuación del tronco `y² + z² = R1²`. Para `α = 90°` coincide con el notebook; verifiqué que con `R1=100, R2=50` da `(0, 50, 86.6)` para `φ=0` y `(±50, 0, 100)` para `φ=±π/2`.

### Paso 1b — Vértices de la curva (`vertices_curva`, `phi_curva_mitad`)

Cuando la derivación está muy inclinada, la curva de intersección **tiene esquinas** (al "desenrollarla" sobre el caño base, la tangente gira de golpe). Estas funciones:

1. Detectan esas esquinas: puntos donde la tangente gira más de `ANGULO_VERTICE` grados. Las refinan en 5 pasadas.
2. Reparten los `n+1` nodos de la curva **equiespaciados en longitud de arco**, pero obligando a que **haya un nodo exactamente en cada vértice**.

Esto es mejor que el reparto uniforme en `φ` del notebook, que no respetaba esquinas y daba elementos de tamaño muy distinto.

### Paso 2 — Malla de la derivación (`malla_derivacion_mitad`)

Es la versión generalizada de la función del notebook (secciones 4.1 a 4.4). Diferencias:

- Los nodos de la curva vienen de `phi_curva_mitad` (con nodo en cada vértice).
- Cada "generatriz" (línea desde la curva hasta la boca) tiene largo `H_2 − t_curva`, donde `t_curva = nodos_curva · eje_s` es la posición axial. Si alguna da ≤ 0 se avisa que `H_2` es demasiado corto.
- La cantidad de elementos `M` a lo largo se elige para que sean ~cuadrados: `M = round(mediana(largo) / tam_s)`.
- Además de crecer geométricamente (`razon`), las generatrices **giran gradualmente** de `φ_curva` a `φ_boca` (la boca es equiespaciada), así la boca queda con nodos parejos aunque la curva no lo sea.

Devuelve los mismos tipos de grupos de nodos: `curva`, `boca`, `simetria`.

### Paso 3 — La corona (`malla_base_borde_mitad`)

Acá aparece algo nuevo respecto del notebook. Una malla de derivación que termina de golpe sobre el caño directo da elementos muy deformados alrededor de la curva. Para evitarlo se agregan `capas_base` (3) capas de cuadriláteros **sobre el caño directo, alrededor de la curva**: la corona (grupo 1).

Cómo se construye:

1. Se **desenrolla** el cilindro base: cada punto pasa a ser `(x, w)` con `w = R1·atan2(z, y)`. En ese plano la curva es una línea y las capas se arman como desplazamientos hacia afuera.
2. Se calcula, para cada nodo de la curva, la **normal hacia afuera**.
3. Cada capa se desplaza `rel × largo del elemento` hacia afuera, para mantener la relación de aspecto.
4. En las **esquinas convexas** (vértices) se agrega un **abanico**: varias direcciones y un bloque de `capas × capas` celdas por sector, para no deformar la malla en el codo de la curva.
5. Se vuelve a 3D y se orientan las normales hacia afuera.

Devuelve, entre otros, `borde_orden`: la cadena ordenada de nodos del borde exterior de la corona, que se usa en el paso 5.

### Paso 4 — Unión (`unir_mallas`)

Une la derivación y la corona **identificando los nodos de la curva** (comparten esos nodos). Verifica que coincidan; si no, error `Los nodos a identificar no coinciden`. Devuelve un `mapa` para renumerar el borde.

### Paso 5 — Completar el caño directo (`completar_base_mitad`)

Falta mallar el resto del caño directo, **desde las bocas hasta la corona**. Se hace sobre el desarrollo del cilindro, dividiéndolo en **tres bloques estructurados**: izquierdo, inferior y derecho (el superior ya lo ocupa la corona).

- Los nodos de las bocas se reparten equiespaciados (`N_directo/2` en la mitad).
- Cada bloque se rellena por **interpolación transfinita** (Coons): los nodos interiores se calculan a partir de los cuatro bordes del bloque. Es la fórmula larga de `p = ...` en el código.
- Se prueban varios valores de `n` (cuántas filas del borde de la corona van a los bloques laterales) y se **elige el que da mejor calidad**.
- Después un **suavizado por optimización** mueve los nodos interiores mientras mejore el peor elemento (las bocas, los planos de simetría y la corona quedan fijos).
- Todo se pasa a 3D y se orientan las normales hacia afuera.

Antes de esto, `generar_malla_tee` revisa que las bocas no corten la corona. Si `H_1` o `H_3` son muy cortos, **los amplía automáticamente** y avisa (`AVISO: la corona llega a x = ...; se amplía H_1 ...`).

### Paso 6 — La cuña (`cuna_simetria`)

En los extremos de la corona, donde toca el plano de simetría, puede quedar un elemento muy grande comparado con el tamaño objetivo. Si supera `umbral_cuna` (1.4) veces el tamaño objetivo, se agrega una **fila de elementos delgados pegada a la línea de simetría** (la "cuña"). Genera algunos **triángulos**, que se guardan como cuadriláteros con el último nodo repetido `[a, b, c, c]`. Al espejar se completan como **rombos**, así que la pieza final no tiene triángulos. Para acomodar los nodos alrededor usa `suavizar_zona`.

### Paso 7 — Alinear (`alinear_malla`)

Endereza las líneas del caño directo y lleva los elementos hacia la forma cuadrada, **sin cambiar la conectividad** (solo mueve nodos). Cada elemento "propone" una posición para cada uno de sus nodos (la que lo haría un paralelogramo perfecto, regla `b + d − c`) y el nodo se mueve hacia el promedio de las propuestas. Un movimiento solo se acepta si la calidad de sus elementos no baja de `piso` (0.5), de modo que **la malla no puede invertirse**. Los nodos del borde, los compartidos entre grupos y los de la cuña quedan fijos.

### Paso 8 — Espejar (`espejar_malla`)

Igual a la del notebook (sección 5) con un detalle más: los **triángulos de la cuña** (con dos nodos sobre el plano de simetría) no se duplican; junto con su espejo forman un rombo `[b, c, a, c']`.

## 23. Control de calidad (`control_malla`)

Devuelve un diccionario con:

| Clave | Qué mide |
|---|---|
| `forma_media`, `forma_min` | Calidad de forma de los elementos (1 = cuadrado perfecto). |
| `libres` | Aristas con un solo elemento (**bordes abiertos**). |
| `mas_de_2` | Aristas compartidas por más de 2 elementos (error). |
| `orientacion_inconsistente` | Elementos vecinos con orientación contraria. |
| `normales_invertidas` | Elementos cuya normal apunta hacia adentro. |

En una pieza sana, `mas_de_2`, `orientacion_inconsistente` y `normales_invertidas` valen 0. `libres` son los bordes de las bocas.

**Resultado real** (ejemplo del docstring: `OD_1=200, OD_2=100, H_1=350, H_2=350, H_3=400, Alpha_2=30`):

| | Valor |
|---|---|
| Media pieza | 3747 nodos, 3604 elementos |
| Pieza completa | **7305 nodos, 7206 elementos** (800 derivación + 130 corona + 6276 resto) |
| `M` (capas en la derivación) | 20 |
| `forma_media` / `forma_min` | 0.907 / 0.337 |
| `libres` | 200 = 80 (boca de ingreso) + 80 (boca de salida) + 40 (boca de derivación) |
| Orientación / normales | 0 inconsistentes, 0 invertidas |
| Tiempo | ~8 s |

Con `Alpha_2 = 90` (derivación perpendicular): 8639 nodos, 8540 elementos, `forma_media` 0.983, `forma_min` 0.666 y ~0.4 s. Es más rápido y mejor porque la curva de intersección es simétrica y sin esquinas.

> Al correr con `Alpha_2=30` aparece un `RuntimeWarning: invalid value encountered in divide` en la línea 445 (cálculo de calidad de elementos del suavizado). Es una división por cero con elementos degenerados (los triángulos); no impide generar la malla, pero conviene saberlo.

## 24. Salida de `generar_malla_tee`

Un diccionario con:

- `nodos`, `elementos`, `grupo`: la pieza completa (ya espejada).
- `divisiones`: `N_directo`, `N_s`, `tam`, `tam_s`, `H_1/H_2/H_3` (**los efectivamente usados**, por si se ampliaron), `M_derivacion`, `capas_base`, etc.
- `control`: resultado de `control_malla`.
- `media_pieza`: nodos, elementos y grupo antes de espejar.
- `ejes`: `{derivacion: eje_s, directo: eje_x}`, útiles para graficar y para `control_malla`.

## 25. Funciones que existen pero no se usan en el flujo principal

- `uniformizar_tamano`: elimina elementos donde la malla está comprimida (fusionando nodos opuestos de un cuadrilátero). Es **opcional**: ni `tee.py` ni `ensamblador.py` la llaman. Los parámetros `tam_dir` y `tam_der` indican el tamaño objetivo.
- `distribucion_geometrica`: es copia de la del notebook y del manual (sección 3).

Dato útil: `main.py` importa `tee` solo para usar `tee.control_malla` y verificar la malla fusionada final de toda la red.

---

---

# Parte 6 — El ensamblador (`src/ensamblador.py`)

Hasta acá tenemos piezas sueltas (Partes 2, 3 y 5) y la geometría de la red (Parte 4). El ensamblador las junta: **genera la malla de cada tramo y accesorio, la orienta y la ubica en el espacio, y después las funde en una sola malla** con nodos compartidos en cada unión.

Tiene dos funciones principales: `ensamblar` (genera las piezas) y `fusionar_piezas` (las une).

## 26. Visión general

```
leer_mermaid ─► calcular_posiciones ─► ensamblar ─► fusionar_piezas
 (nodos, tramos)   (pos, dir_in, ref_in)   (lista de piezas)   (malla única + soldaduras)
```

`ensamblar(nodos, tramos, N)` devuelve una **lista de piezas**, cada una `{"nodos", "elementos", "nombre"}`. Hay una pieza por cada tramo recto y una por cada accesorio (tee, codo, reducción). Los nodos `Libre` no generan pieza.

Resultado real con `Caso_1.md` y `N=48`:

| Pieza | Nodos | Elementos |
|---|---|---|
| tramo 1→2 | 3936 | 3888 |
| tee 2 | 2206 | 2143 |
| tramo 2→3 | 8352 | 8304 |
| tramo 2→6 (derivación) | 2624 | 2592 |
| codo 3 | 864 | 816 |
| tramo 3→4 | 4464 | 4416 |
| reducción 4 | 1488 | 1440 |
| tramo 4→5 | 5568 | 5520 |

## 27. El recorrido y cómo se ubica cada accesorio

`ensamblar` recorre el árbol igual que `red.py` (pila, desde la raíz). Para cada nodo:

1. Toma `pos`, `dir_in` y `ref_in` de `calcular_posiciones`.
2. Si es un accesorio, genera su pieza con la función correspondiente (`_pieza_codo`, `_pieza_reduccion`, `_pieza_tee`).
3. Por cada tramo que sale del nodo, genera el tubo (`_pieza_tramo`).

**Dónde se coloca cada accesorio** (esto es lo que hay que entender bien, porque `pos` es el nodo del grafo):

| Accesorio | Punto de inserción |
|---|---|
| **Reducción** | `pos` es su punto de **entrada**. La pieza se extiende `L` hacia adelante. |
| **Tee** | `pos` es el punto donde **se cruzan** tronco y derivación (el origen local `x = 0` de `tee.py`). El tronco se extiende `H_1` hacia atrás y `H_3` hacia adelante. |
| **Codo** | `pos` es el **vértice**: la intersección de las rectas de entrada y de salida, prolongadas. La malla del codo en cambio arranca en el punto de tangencia, retrocedido `T` desde el vértice. |

La distancia `T` es la geometría estándar de un codo (función `_tangente_codo`):

```
T = R · tan(Angulo / 2)
```

Para `R = 225` y `Angulo = 45°`: `T = 93.2`. Es lo que se retrocede para ubicar el inicio del codo, y también lo que va a **cubrir** el codo en cada extremo (sección 29).

## 28. Orientar cada pieza con `alinear_pieza`

Como vimos en la sección 11.3, `alinear_pieza` rota una pieza para que su entrada y su eje normal calcen con los globales. Acá se usan:

```python
EJE_ENTRADA_LOCAL = (1, 0, 0)    # las 4 piezas entran por +X
EJE_NORMAL_LOCAL  = (0, -1, 0)   # eje normal del desvío, igual para codo y tee
```

El eje normal es **−Y y no +Y**. El comentario del código explica por qué: `codo.py` y `tee.py` curvan su tangente de `(1, 0, 0)` hacia `(cos t, 0, sin t)`, y esa trayectoria corresponde a una rotación de Rodrigues alrededor de **−Y**; se verificó comparando la tangente real de salida de la malla con la dirección calculada por `red.py` (con +Y salía invertida).

Qué eje **global** se pasa a `alinear_pieza` según la pieza:

| Pieza | `eje_entrada_global` | `eje_normal_global` |
|---|---|---|
| Codo y Tee | `dir_in` | `_eje_desvio(ref_in, dir_in, Giro)` (el mismo que usa `red.py`) |
| Reducción | `dir_in` | `ref_in` |
| Tramo recto | `dir_tramo` | `ref_tramo` |

Después de rotar se **traslada** sumando la posición (`nodos_rot + pos`).

La Tee y el Codo usan el mismo `_eje_desvio` que `red.py`: así la dirección en la que **realmente** sale la malla coincide con la dirección en la que `red.py` ubicó el nodo siguiente.

## 29. Los tramos rectos y sus "recortes"

Esta es la respuesta a la pregunta que quedó abierta en la Parte 4. En el grafo, `L` es la distancia entre centros de nodos, pero **los accesorios ya ocupan parte de esa distancia** con su propia malla. Si el tubo midiera `L` completo, se superpondría con ellos. Por eso:

```python
L_tubo = L − recorte_inicio − recorte_fin
pos_tubo = pos_origen + dir_tramo · recorte_inicio
```

El tubo **arranca más adelante** (donde termina el accesorio de origen) y **termina antes** (donde empieza el accesorio de destino).

Cuánto cubre cada accesorio:

| Accesorio | Recorte en su salida (`recorte_inicio`) | Recorte en su entrada (`recorte_fin`) |
|---|---|---|
| Tee, tramo `S` (directo) | `H_3` | `H_1` |
| Tee, tramo `D` (derivación) | `H_2` | — |
| Codo | `T = R·tan(Angulo/2)` | `T` |
| Reducción | `L` de la reducción | 0 |
| Libre | 0 | 0 |

Ejemplo con `Caso_1` (largos de tubo resultantes):

| Tramo | `L` | Recorte inicio | Recorte fin | Tubo real |
|---|---|---|---|---|
| 1→2 | 1000 | 0 (Libre) | 200 (`H_1` tee) | **800** |
| 2→3 | 2000 | 211 (`H_3` tee) | 93.2 (`T` codo) | **1695.8** |
| 2→6 | 1000 | 200 (`H_2` tee) | 0 (Libre) | **800** |
| 3→4 | 1000 | 93.2 (`T` codo) | 0 (reducción) | **906.8** |
| 4→5 | 1000 | 250 (`L` reducción) | 0 (Libre) | **750** |

Si los recortes suman más que `L` el programa para con un error que pide revisar `L`, `H_1/H_2/H_3`, `R`, `Angulo` o el largo de la reducción.

## 30. El número de divisiones `N` en la derivación

Toda la red usa el `N` que se pasa al programa (`N=48` en el ejemplo)... **excepto la derivación de una tee**. La tee decide su propio `N_s` (el que da elementos casi cuadrados en la derivación, ver sección 20) y **ese es el `N` que se usa de ahí en adelante en toda la rama**: el tramo siguiente y cualquier accesorio posterior.

Para eso el ensamblador lleva dos diccionarios:
- `N_in[nodo]`: el `N` con el que se genera lo que **entra** a cada nodo.
- `OD_in[nodo]`: el diámetro del tramo que entra a cada nodo. Codo y Tee no tienen `OD` en sus propios parámetros (en el Mermaid solo figura `OD2` para la tee), así que **toman el diámetro del tramo que les llega** (`OD_1` de la tee, `OD` del codo).

En `Caso_1` con `N=48`, el programa imprime: `Tee 2: tronco N=48 -> derivacion N_s=32`. El tramo 2→6 se genera entonces con N=32.

## 31. Que los nodos coincidan: fase y verificaciones

Para FEM no alcanza con que dos piezas **se vean** juntas: los nodos de una y otra deben **coincidir exactamente**. El código toma tres precauciones, cada una con un error explícito si no se cumple.

### 31.1 Fase de los tubos

Un tubo es simétrico en geometría, pero su malla discreta no: los `N` nodos del perímetro están en ángulos fijos. Para que el anillo de un tubo caiga sobre el del accesorio vecino, hay que girar el tubo a la **misma fase** (el mismo `ref`) que el accesorio. Por eso `_pieza_tramo` recibe `ref_tramo`:

- Si el destino es Codo o Tee: `ref_tramo = _eje_desvio(ref_in, dir_in, Giro)` del destino (esas piezas son un solo cuerpo rígido, su cara de entrada ya está girada `Giro`).
- Si es Reducción o Libre: `ref_tramo = ref_in` sin modificar.

### 31.2 H_1 y H_3 de la tee tienen que ser suficientes

Recordá que `tee.py` **amplía en silencio** `H_1` y `H_3` si la corona no entra (sección 22, paso 5). Pero los tubos vecinos se recortaron con los valores que escribió el usuario. Resultado: la tee terminaría en otro lugar que el tubo.

Por eso `_pieza_tee` compara los `H_1/H_3` reales con los nominales y, si difieren, **lanza un error** que dice qué valor mínimo declarar: `hace falta declarar H1>=... y H3>=...`. Esto explica el `H3=211` del `Caso_1`: es el mínimo que necesita esa tee (con `Alpha=60`).

### 31.3 Dos accesorios enchufados directo con un tramo en el medio

Si un Codo/Tee conecta con otro Codo/Tee por un tramo recto, ese tramo tiene una **única fase** en toda su longitud, pero necesita calzar con las dos caras. Eso solo es posible si la diferencia de fase entre ellas es **múltiplo del paso angular** `360°/N`. Si no, el programa avisa cuánto desfase hay y propone tres salidas: cambiar el `Giro` de uno de los accesorios, cambiar `N`, o intercalar una reducción o un nodo libre.

## 32. Fusión de piezas (`fusionar_piezas`)

Convierte la lista de piezas en **una sola malla** con nodos compartidos.

1. **Concatenar.** Se apilan todos los nodos y se suman *offsets* a los índices de los elementos.
2. **Detectar nodos "iguales".** Se redondean las coordenadas a la tolerancia `tol` (1e-3):

```python
clave = np.round(nodos / tol).astype(int)
_, idx_unico, inverso = np.unique(clave, axis=0, return_index=True, return_inverse=True)
```

   `np.unique` con `axis=0` le asigna el mismo id a todos los nodos con la misma clave.
3. **Reasignar.** `nodos_fusionados = nodos[idx_unico]` y `elementos_fusionados = inverso[elementos]`.

> Una limitación a tener en cuenta: el redondeo agrupa por celdas. Dos nodos muy cercanos pero separados por el borde de una celda **podrían no fusionarse**. Con `tol = 1e-3` y coordenadas del orden de cientos de mm no es un problema práctico, pero es la razón de que la coincidencia deba ser casi exacta (por eso toda la preocupación por la fase).

Resultado en `Caso_1`: **29 502 nodos sueltos → 29 182 fusionados** (se eliminaron 320) y **29 119 elementos**.

## 33. Zonas de soldadura

Un nodo fusionado que recibió **2 o más nodos originales** es un nodo donde dos piezas se unen: un **nodo de soldadura**.

Pero un anillo de unión tiene varios nodos, y hay que saber cuáles pertenecen al **mismo anillo**. `_agrupar_soldadura` lo resuelve con **Union-Find** (conjuntos disjuntos): recorre las aristas de todos los elementos y, si los dos extremos son nodos de soldadura, los une en la misma zona. Cada zona queda como un array ordenado de ids; las zonas se ordenan por su primer id para que la salida sea reproducible.

En `Caso_1` salen **7 zonas**, correspondientes a las 7 uniones entre piezas:

| Zona | Unión | Nodos |
|---|---|---|
| 1 | tramo 1→2 ↔ tee (ingreso) | 48 |
| 2 | tee ↔ tramo 2→6 (derivación) | 32 |
| 3 | tee ↔ tramo 2→3 (salida directa) | 48 |
| 4 | tramo 2→3 ↔ codo | 48 |
| 5 | codo ↔ tramo 3→4 | 48 |
| 6 | tramo 3→4 ↔ reducción | 48 |
| 7 | reducción ↔ tramo 4→5 | 48 |

Suma: 320 nodos, igual a los que se eliminaron al fusionar. ✔

(El orden de la tabla es el de la lista de uniones; en el programa las zonas salen ordenadas por id.)

## 34. Verificación de la malla final

`main.py` llama a `tee.control_malla` sobre la malla fusionada. Para `Caso_1` con `N=48`:

```
forma_media = 0.994, forma_min = 0.504
libres = 128, mas_de_2 = 0, orientacion_inconsistente = 0
```

- `libres = 128` son los bordes abiertos que **deben** quedar: la boca inicial (48), el extremo de la derivación (32) y el extremo final (48).
- `mas_de_2 = 0` y `orientacion_inconsistente = 0` indican que las uniones quedaron bien soldadas y que todas las normales apuntan al mismo lado. Es la comprobación de que el trabajo de fase y recortes funcionó.

---

---

# Parte 7 — Exportar, graficar y la interfaz gráfica

Con las partes anteriores ya tenemos la malla única fusionada. Esta parte explica cómo se **guarda** (`exportador.py`), cómo se **dibuja** (`graficador.py`) y cómo se **usa todo desde una ventana** (`mesh_preforma_gui.py`).

## 35. El exportador (`src/exportador.py`)

Guarda la malla en un archivo de texto plano `.malla` que otro programa (el de FEM) puede leer. Tiene dos funciones: `escribir_malla` y `leer_malla`.

### 35.1 Formato del archivo

Cuatro secciones, cada una con un encabezado `### ... ###`:

```
### MERMAID ###
# ```mermaid
# flowchart TD
# 1((Libre))
# 2["Tee:OD2=100,Alpha=60,..."]
# ...
### NODOS ###
N 7291
0 0.000000 -75.000000 0.000000
1 0.000000 -72.444437 -19.411428
...
### ELEMENTOS ###
E 7260
15 17 41 39
17 19 43 41
...
### SOLDADURAS ###
Z 7
0 24 984 985 986 ... 
1 16 1301 1306 1307 ...
```

| Sección | Contenido |
|---|---|
| `MERMAID` | El diagrama de origen, **cada línea comentada con `# `**. Así el archivo guarda con qué red se generó. |
| `NODOS` | `N <cantidad>`, y después una línea por nodo: `id x y z` (6 decimales). |
| `ELEMENTOS` | `E <cantidad>`, y después una línea por elemento con 4 ids de nodo. |
| `SOLDADURAS` | `Z <cantidad de zonas>`, y una línea por zona: `número cantidad id1 id2 ...`. |

Detalles importantes:
- **Los ids de nodo empiezan en 0** y son consecutivos.
- Un **triángulo** se guarda con el último nodo repetido: `[a, b, c, c]`. Es la misma convención de todo el proyecto.
- La sección `SOLDADURAS` es **opcional**: si se llama `escribir_malla` sin `zonas_soldadura`, no se escribe.
- Cada zona de soldadura es un anillo (una unión entre piezas), tal como las calculó `fusionar_piezas` (sección 33).

### 35.2 `escribir_malla(path, nodos, elementos, mermaid_texto, zonas_soldadura)`

Convierte los arreglos a `float` e `int` y escribe en el orden de la tabla. No tiene lógica: es solo formateo.

### 35.3 `leer_malla(path)`

Hace el camino inverso y devuelve `(nodos, elementos, mermaid_texto, zonas_soldadura)`:
- Verifica los encabezados (con `assert`).
- Quita el `# ` de cada línea del Mermaid para recuperar el texto original.
- Lee las zonas solo si existe la sección `SOLDADURAS`; si no (formato viejo), devuelve `[]`.

Una particularidad: al leer los nodos usa la **posición de la línea** y no el número de id (ignora la primera columna). Por eso los ids tienen que ir en orden, empezando por 0, que es como los escribe `escribir_malla`.

### 35.4 Verificación

Probé el ciclo completo con `Caso_1` y `N=24`: escribir y volver a leer.

| Prueba | Resultado |
|---|---|
| Tamaño del archivo | ~420 KB (7291 nodos, 7260 elementos, 7 zonas) |
| Diferencia máxima en coordenadas | 5·10⁻⁷ (por los 6 decimales) |
| Elementos | idénticos |
| Zonas de soldadura | 7, idénticas |
| Texto del Mermaid | recuperado igual |
| Archivo sin soldaduras | se lee bien, devuelve `[]` |

Dato curioso: con `N=24` la derivación de esa tee usa `N_s = 16`, y por eso su zona de soldadura tiene 16 nodos (en vez de los 32 que vimos con `N=48`).

## 36. El graficador (`src/graficador.py`)

Dibuja con matplotlib. **No sabe nada de cañerías**: solo recibe arreglos de nodos `(N, 3)` y elementos `(E, 4)`. Tiene cuatro funciones.

### 36.1 `graficar_malla`

Dibuja la malla en 3D como una superficie de cuadriláteros (`Poly3DCollection`).

- `grupo` (opcional): colorea cada elemento según su grupo. Los colores por defecto son los de la tee: 0 = derivación (naranja), 1 = corona (verde), 2 = directo (celeste).
- `colores`: diccionario para cambiar esos colores. En `main.py` se pasa `{0: "#ffcf80"}`, un único color para toda la malla.
- `ax` y `mostrar`: permiten dibujar varias cosas sobre los mismos ejes antes de mostrar (`mostrar=False`).
- Ajusta los límites de los ejes con **la misma escala en X, Y y Z** (cubo `set_box_aspect((1,1,1))`), para que no se deforme la geometría.

### 36.2 `graficar_desarrollado`

Dibuja la malla **"desenrollada" en 2D**: `(x, R·atan2(z, y))`. Es el mismo desarrollo que se usa dentro de `tee.py` (sección 22) y sirve para ver de un vistazo la regularidad de los elementos sobre un caño. Si hay una derivación en otro eje, esos elementos salen deformados en este gráfico.

### 36.3 `graficar_soldaduras`

Marca sobre la malla las zonas de soldadura como **anillos de color** (magenta por defecto):

- Los nodos de una zona llegan **ordenados por id**, no por posición en la vuelta. Para dibujarlos como una línea cerrada hay que ordenarlos por ángulo: lo hace `_ordenar_anillo`.
- `_ordenar_anillo` encuentra el plano del anillo con una **descomposición SVD** (la dirección de menor varianza es la normal), arma una base 2D dentro de ese plano y ordena con `arctan2`.
- El anillo se dibuja un 8 % **agrandado** (`expandir=1.08`) para que quede como relieve por fuera del caño. Esto es porque matplotlib no hace un z-buffer real en 3D y, si no, la malla lo taparía.

### 36.4 `graficar_ejes`

Dibuja el **esqueleto de la red**: un segmento por cada tramo y un punto en cada nodo, con etiqueta `id:Tipo` y color según el tipo (Libre gris, Tee naranja, Codo azul, Reducción verde). Se usa para verificar el ruteo de `red.py` **antes** de mallar, y en `main.py` se superpone sobre la malla.

## 37. La interfaz gráfica (`mesh_preforma_gui.py`)

Es una ventana de escritorio hecha con **tkinter**, con matplotlib incrustado. Se ejecuta con:

```
uv run python mesh_preforma_gui.py [archivo.md]
```

Si no se pasa archivo, arranca con una **plantilla de ejemplo** (tee + codo + reducción).

> Nota: no pude ejecutar la interfaz en mi entorno (no tiene entorno gráfico ni tkinter). Esta sección sale de leer el código, que sí verifiqué que compila. Todas las funciones de fondo (`ensamblar`, `fusionar_piezas`, `escribir_malla`) son las mismas que ya ejecuté en las partes anteriores.

### 37.1 Qué se ve

- **Izquierda:** un editor de texto con el Mermaid.
- **Derecha:** dos pestañas, *Diagrama Mermaid* (cajas y flechas en 2D) y *Vista 3D* (traza o malla), más un panel **Estado** (log) abajo.
- **Barra superior:** Nuevo archivo, Abrir, Guardar Mermaid, Deshacer/Rehacer, el spinbox de **N** (por defecto 24, rango 4 a 400), **Previsualizar malla**, **Mallar** y **Guardar .malla**.
- **Segunda barra:** el combo *Agregar desde nodo* y los botones **+ Tramo recto, + Tee, + Codo, + Reducción**.

### 37.2 Previsualizar y Mallar

Ambos botones llaman a la misma función, `_generar_malla`, que hace exactamente lo que hace `main.py`:

1. Parsea el texto del editor (`parsear_mermaid`).
2. `ensamblar` (con el `N` del spinbox) y `calcular_posiciones`.
3. `fusionar_piezas` y `tee.control_malla`.
4. Dibuja malla + soldaduras + ejes en la pestaña 3D.
5. Escribe en el log el control de malla y la cantidad de zonas de soldadura.

La única diferencia: **Mallar** muestra una advertencia emergente si `mas_de_2` u `orientacion_inconsistente` no valen 0; **Previsualizar** no, para no interrumpir mientras se está iterando.

Si algo falla (por ejemplo los errores de `H1/H3` de la sección 31), se muestra un cuadro de error con el mensaje y se escribe en el log.

Los avisos que el programa imprime con `print` (por ejemplo `Tee 2: tronco N=48 -> derivacion N_s=32`) se **capturan** (`redirect_stdout`) y se muestran en el panel Estado.

### 37.3 Refresco automático

Mientras se escribe, la ventana se actualiza sola, sin apretar nada:

- Cada cambio en el editor dispara `_on_editor_modificado`, que **invalida cualquier mallado anterior** (así no se puede guardar un `.malla` que ya no corresponde al texto) y programa un refresco.
- El refresco usa **debounce**: se cancela y reprograma con cada tecla, y recién corre cuando pasan **500 ms sin cambios**.
- Redibuja el diagrama 2D y una **traza liviana** en 3D (solo los ejes, sin mallar), que es rápida.
- Si el texto tiene un error de sintaxis, no muestra mensajes: simplemente no actualiza. Los errores se muestran recién al apretar Previsualizar o Mallar.

Como consecuencia, **después de editar hay que volver a apretar Mallar** antes de poder guardar.

### 37.4 Guardar

- **Guardar .malla…**: si todavía no se generó la malla (botón Mallar), avisa. Si no, abre un diálogo (nombre sugerido: el del `.md` cargado con extensión `.malla`) y llama a `escribir_malla` con el texto actual del editor como encabezado y las zonas de soldadura.
- **Guardar Mermaid…**: guarda el texto del editor como `.md`.

### 37.5 Los botones de ayuda para escribir la red

Son el aporte más práctico de la interfaz: arman el Mermaid sin tener que escribir la sintaxis a mano. Cada uno abre un formulario (`_pedir_valores`) con los datos numéricos y **inserta las líneas** en el lugar correcto del texto (`_insertar_lineas`).

Para decidir qué hacer, necesitan saber el **grado** del nodo de origen (cuántos tramos lo tocan), porque la topología exige grados fijos (Libre 1, Codo 2, Reducción 2, Tee 3; sección 16):

| Situación del nodo de origen | Qué hace el botón |
|---|---|
| **Grado 0** (nodo suelto, el primero de la red) | Crea el accesorio, un **tramo de entrada** que lo alimenta y los nodos libres de salida. Por eso el formulario pide además largo y diámetro del tramo de entrada. |
| **Grado ≥ 1 con tramo entrante** (un extremo libre al final de la red) | **Convierte** ese nodo `Libre` en el accesorio (reemplaza su línea de declaración) y agrega los tramos de salida. Así no queda un `Libre` con 2 tramos, que `red.py` rechazaría. |
| **Grado ≥ 1 sin tramo entrante** | Avisa que no es un extremo disponible. |

Por tipo:
- **+ Tramo recto**: solo sirve para **iniciar** la red desde un nodo suelto. Si el nodo ya tiene un tramo, avisa que para alargar hay que editar el `L` directamente.
- **+ Tee**: pide `OD2, Alpha, H1, H2, H3, Giro` y los largos y diámetros de los tramos `S` (directo) y `D` (derivación).
- **+ Codo**: pide `R, Angulo, Giro` y largo y diámetro del tramo de salida.
- **+ Reducción**: pide `OD1, L, OD2` y el largo del tramo de salida. El tramo de entrada toma el diámetro `OD1` y el de salida `OD2`, porque la reducción cambia justo de uno al otro.

El diámetro de los formularios se **sugiere** a partir del último tramo escrito (`_od_sugerido`, 150 si no hay ninguno), y después de agregar algo, el combo de origen pasa al nuevo extremo, para poder encadenar piezas rápido.

### 37.6 Deshacer y diagrama

- Deshacer/Rehacer (Ctrl+Z / Ctrl+Y) usan el historial del editor. Los reemplazos de texto completos (como los botones `+`) cuentan como **un solo paso** de deshacer.
- El diagrama 2D se dibuja con matplotlib (cajas con `FancyBboxPatch`, flechas con `annotate`). `_layout_mermaid` ubica los nodos en niveles por profundidad desde la raíz, como un `flowchart TD`.

## 38. Mapa final: de punta a punta

```
Editor / archivo .md                (texto Mermaid)
   │ mermaid_reader.parsear_mermaid
   ▼
nodos + tramos                      (grafo)
   │ red.calcular_posiciones
   ▼
pos, dir_in, ref_in por nodo        (geometría global)
   │ ensamblador.ensamblar          ── tubo.py / codo.py / reduccion.py / tee.py
   ▼
lista de piezas                     (una malla por tramo y accesorio, ya orientadas)
   │ ensamblador.fusionar_piezas
   ▼
malla única + zonas de soldadura
   │ tee.control_malla              → verificación
   │ exportador.escribir_malla      → archivo .malla
   │ graficador.*                   → vista 3D
```

## 39. Comentarios que estaban desactualizados (ya corregidos)

Al leer el código encontré comentarios que no coincidían con lo que el programa hace. No afectaban el funcionamiento y **ya los corregí**:

| Archivo | Decía | Ahora dice |
|---|---|---|
| `mesh_preforma_gui.py` (docstring) | `uv run python app_gui.py` | `uv run python mesh_preforma_gui.py` |
| `ensamblador.py` (docstring del módulo) | Que la fusión de nodos era un paso pendiente | Describe el flujo real: `ensamblar()` genera las piezas y `fusionar_piezas()` las une y devuelve las zonas de soldadura |
| `main.py` (comentario) y `README.md` | Soldaduras "en rojo" | "en magenta", que es el color por defecto de `graficar_soldaduras` |
| `main.py` (docstring) | Mencionaba `prueba_red.md` | `ejemplos/Caso_1.md`, el único ejemplo que existe en esa carpeta |

Antes ya se habían corregido los docstrings de `codo.py` (la entrada es por +X) y de `reduccion.py` (x=0 y x=L).

---

## Fin del recorrido

Ya están explicadas todas las piezas del programa: tubo, codo, reducción y tee (Partes 2, 3 y 5), la lectura del Mermaid y la geometría de la red (Partes 3 y 4), el ensamblado y la fusión (Parte 6), y la exportación, los gráficos y la interfaz (Parte 7).

Un buen siguiente paso para practicar: cargar `ejemplos/Caso_1.md` en la interfaz, cambiar `Giro`, `Alpha` y `N`, y mirar cómo cambia la malla y qué errores de las secciones 31 y 22 aparecen.

---

*Este manual se va ampliando a medida que se agregan pasos al programa.*
