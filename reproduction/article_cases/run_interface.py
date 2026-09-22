#!/usr/bin/env python3
"""Reproduce one- and four-benzene graphene interface parameter datasets."""
from pathlib import Path
import json
from ase.io import read
from hopping3d_eht.interface import parameterize_interface
from hopping3d_eht.multi_interface import parameterize_coverage_interface

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent/'results'; OUT.mkdir(exist_ok=True)
one=parameterize_interface(read(ROOT/'examples/article_systems/graphene_benzene/graphene_benzene.cif'))
four=parameterize_coverage_interface(read(ROOT/'examples/article_systems/graphene_benzene/graphene_4benzene.cif'))
(OUT/'graphene_benzene_single.json').write_text(json.dumps(one,indent=2)+'\n')
(OUT/'graphene_benzene_coverage.json').write_text(json.dumps(four,indent=2)+'\n')
assert one['detection']['molecule_type']=='benzene'
assert four['detection']['n_adsorbates']==4
print('Wrote interface parameter files to',OUT)
