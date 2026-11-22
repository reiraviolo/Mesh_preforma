"""Estructuras de datos compartidas por la lectura de entrada, los accesorios y el ensamble."""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Tramo:
    id: int
    dir: np.ndarray  # (3,) vector unitario de Nodo_I hacia Nodo_F
    long: float  # distancia entre nodos (centro a centro de accesorio)
    od: float
    nodo_i: int
    nodo_f: int


@dataclass
class Nodo:
    id: int
    tipo: str | None  # "Tee", "Codo", "Reduccion" o None (extremo libre)
    params: dict[str, float] = field(default_factory=dict)
    pos: np.ndarray | None = None  # se completa en el ensamble


@dataclass
class Puerto:
    centro: np.ndarray  # centro del anillo de la boca
    direccion: np.ndarray  # unitario, hacia afuera del accesorio
    od: float
    ids_anillo: np.ndarray  # IDs de los nodos de malla de esa boca


@dataclass
class Malla:
    nodos: np.ndarray  # (n, 3)
    elementos: np.ndarray  # (m, 4)
    puertos: list[Puerto] = field(default_factory=list)


@dataclass
class Red:
    tramos: dict[int, Tramo]
    nodos: dict[int, Nodo]
    params: dict[str, float] = field(default_factory=dict)  # globales, ej. N

    def tramos_de(self, id_nodo: int) -> list[Tramo]:
        """Tramos que llegan o salen del nodo."""
        return [t for t in self.tramos.values() if id_nodo in (t.nodo_i, t.nodo_f)]
