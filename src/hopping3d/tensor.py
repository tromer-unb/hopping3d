#!/usr/bin/env python3
"""Utilities for analyzing symmetric 3D diffusion/mobility tensors."""
from __future__ import annotations
import numpy as np


def symmetrize_tensor(T):
    A=np.asarray(T,dtype=float).reshape(3,3)
    return 0.5*(A+A.T)


def principal_decomposition(T):
    """Return eigenvalues/eigenvectors sorted fastest -> slowest."""
    A=symmetrize_tensor(T)
    vals, vecs=np.linalg.eigh(A)
    order=np.argsort(vals)[::-1]
    return vals[order], vecs[:,order]


def dual_unit_vectors(cell):
    """Unit normals to the (100),(010),(001) lattice planes.

    ASE cell vectors are rows. Columns of inv(cell) are the reciprocal/dual
    vectors without the conventional 2*pi factor.
    """
    A=np.asarray(cell,dtype=float).reshape(3,3)
    inv=np.linalg.inv(A)
    G=np.column_stack([inv[:,i]/np.linalg.norm(inv[:,i]) for i in range(3)])
    return G


def direct_unit_vectors(cell):
    A=np.asarray(cell,dtype=float).reshape(3,3)
    return np.column_stack([A[i]/np.linalg.norm(A[i]) for i in range(3)])


def directional_value(T,n):
    A=symmetrize_tensor(T); v=np.asarray(n,dtype=float); v=v/np.linalg.norm(v)
    return float(v @ A @ v)


def angle_deg(a,b):
    a=np.asarray(a,dtype=float); b=np.asarray(b,dtype=float)
    a=a/np.linalg.norm(a); b=b/np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(a,b)),0.0,1.0))))
