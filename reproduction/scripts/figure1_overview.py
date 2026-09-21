#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from hopping3d.core import load_structure
from ase.data import atomic_numbers,covalent_radii
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; F=ROOT/'generated'; F.mkdir(exist_ok=True)
C=load_structure(ROOT/'inputs/stacked_graphene_4layers.cif')
BN=load_structure(ROOT/'inputs/BN_bulk.cif',repeat=(2,2,2))
WO=load_structure(ROOT/'inputs/W2O6.cif')
fig=plt.figure(figsize=(10.5,3.4))
ax=fig.add_subplot(131); p=C.positions
ax.scatter(p[:,0],p[:,2],s=12,alpha=.7); ax.set(xlabel='x (Å)',ylabel='z (Å)',title='(a) Controlled layered network'); ax.grid(alpha=.15)
for k,(atoms,title) in enumerate([(BN,'(b) Real BN crystal'),(WO,'(c) Chemistry-aware W$_2$O$_6$')],start=2):
    ax=fig.add_subplot(130+k,projection='3d'); pos=atoms.positions; sy=np.asarray(atoms.get_chemical_symbols())
    for sp in sorted(set(sy.tolist())):
        ids=np.where(sy==sp)[0]; size=max(float(covalent_radii[atomic_numbers[sp]])*25,12)
        ax.scatter(pos[ids,0],pos[ids,1],pos[ids,2],s=size,alpha=.7,label=sp)
    ax.set(xlabel='x',ylabel='y',zlabel='z',title=title); ax.legend(frameon=False,fontsize=7)
fig.suptitle('Dimensionality, tensor transport, and chemical heterogeneity',y=.99)
fig.tight_layout(); fig.savefig(F/'figure1_system_overview.png',dpi=400,bbox_inches='tight'); plt.close(fig)
