import numpy as np
from ase import Atoms
from hopping3d.core import direction_vector, build_finite_graph, apply_random_vacancies

def test_cartesian_directions():
    a=Atoms('H',positions=[[0,0,0]])
    assert np.allclose(direction_vector('X',a),[1,0,0])
    assert np.allclose(direction_vector('Y',a),[0,1,0])
    assert np.allclose(direction_vector('Z',a),[0,0,1])

def test_finite_graph_and_reproducible_vacancies():
    a=Atoms('H4',positions=[[0,0,0],[1,0,0],[3,0,0],[4,0,0]])
    g=build_finite_graph(a,1.1)
    assert [len(x) for x in g]==[1,1,1,1]
    b1=apply_random_vacancies(a,.25,seed=4); b2=apply_random_vacancies(a,.25,seed=4)
    assert np.allclose(b1.positions,b2.positions)
