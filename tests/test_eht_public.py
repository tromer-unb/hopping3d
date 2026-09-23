from pathlib import Path
from ase.io import read
from hopping3d_eht.eht import graphene_pi_eht
from hopping3d_eht.article_materials import bn_pi_eht,w2o6_atom_site_eht
from hopping3d_eht.molecular import benzene_monomer,conjugated_pi_dimer_couplings

ROOT=Path(__file__).resolve().parents[1]

def test_graphene_eht_public():
    from runpy import run_path
    cif=ROOT/'examples/paper/01_carbon_stack/stacked_graphene_4layers.cif'
    if not cif.exists(): run_path(str(ROOT/'examples/paper/01_carbon_stack/build_layered_carbon.py'),run_name='__main__')
    a=read(cif); out=graphene_pi_eht(a,cutoff_A=3.75)
    assert out['Heff'].shape==(len(a),len(a))

def test_bn_eht_public():
    a=read(ROOT/'examples/paper/02_bn_tensor/BN_bulk.cif').repeat((2,2,2)); out=bn_pi_eht(a)
    assert out['median_interlayer_rate_scale']>0

def test_w2o6_public():
    a=read(ROOT/'examples/paper/03_w2o6_chemistry/W2O6.cif'); out=w2o6_atom_site_eht(a)
    assert out['edge_counts']['W-W']==0

def test_benzene_electron_hole_channels():
    a=benzene_monomer(); b=a.copy(); b.positions += [0,0,3.4]; d=a+b
    out=conjugated_pi_dimer_couplings(d,range(12),range(12,24),K=2.417879719037213,inter_cutoff_A=8.0)
    assert out['J_hole_eV']>0 and out['J_electron_eV']>0
