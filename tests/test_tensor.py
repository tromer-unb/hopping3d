import numpy as np
from hopping3d.tensor import principal_decomposition, symmetrize_tensor

def test_principal_decomposition_sorted_and_orthogonal():
    A=np.array([[3.,.2,0.],[.2,2.,0.],[0.,0.,1.]])
    vals,vecs=principal_decomposition(A)
    assert np.all(np.diff(vals) <= 0)
    assert np.allclose(vecs.T@vecs,np.eye(3),atol=1e-12)
    assert np.allclose(symmetrize_tensor(A),A)
