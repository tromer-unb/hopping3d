"""Trajectory-resolved kinetic Monte Carlo for 2D/3D atomic hopping networks."""
from .core import *  # noqa: F401,F403
from .graph import active_adjacency, reachable_from_targets, source_drain_reachability
from .survival import simulate_detailed, arrival_curve, adaptive_first_passage
from .tensor import principal_decomposition, directional_value, dual_unit_vectors, symmetrize_tensor

__version__ = "2.1.0"
