#!/usr/bin/env python3
"""Build the controlled four-layer C192 network used in the dimensional-rescue ensemble."""
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import write

A, B, C = 17.04337996, 7.38, 16.4
x_base = np.array([1/24, 1/6, 5/24, 1/12], float)
y_pairs = [(1/6, 0.0), (1/2, 1/3), (5/6, 2/3)]
positions = []
for layer in range(4):
    xshift = (0.71/A) if layer % 2 else 0.0
    z = 3.35 * layer
    for block in range(4):
        for y_top, y_bottom in y_pairs:
            ys = [y_top, y_bottom, y_top, y_bottom]
            for x, y in zip(x_base + 0.25*block + xshift, ys):
                positions.append([float((x % 1.0)*A), float(y*B), z])
atoms = Atoms('C192', positions=positions, cell=[A,B,C], pbc=True)
out = Path(__file__).with_name('stacked_graphene_4layers.cif')
write(out, atoms)
print(out)
