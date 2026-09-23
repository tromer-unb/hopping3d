#!/usr/bin/env python3
"""Graph-theoretic diagnostics used to separate connectivity from KMC kinetics."""
from __future__ import annotations
import numpy as np
from .core import direction_vector, electrode_indices, unit_vector


def active_adjacency(graph, layer_normal, eta, intralayer_cutoff_A=1.9,
                     interlayer_cutoff_A=3.8, interlayer_threshold_A=1.0):
    n = len(graph); adj=[set() for _ in range(n)]; normal=unit_vector(layer_normal)
    for i, edges in enumerate(graph):
        for e in edges:
            inter=abs(float(np.dot(e.dr,normal))) > float(interlayer_threshold_A)
            if inter:
                active=(float(eta)>0.0 and e.r<=float(interlayer_cutoff_A))
            else:
                active=(e.r<=float(intralayer_cutoff_A))
            if active:
                j=int(e.j); adj[i].add(j); adj[j].add(i)
    return adj


def reachable_from_targets(adj, targets):
    reached=set(int(x) for x in targets); stack=list(reached)
    while stack:
        i=stack.pop()
        for j in adj[i]:
            if j not in reached:
                reached.add(j); stack.append(j)
    return reached


def source_drain_reachability(atoms, graph, direction='A', electrode_width_fraction=0.08,
                              layer_normal=None, eta=0.0, intralayer_cutoff_A=1.9,
                              interlayer_cutoff_A=3.8, interlayer_threshold_A=1.0):
    d=direction_vector(direction,atoms); u=np.asarray(atoms.positions)@d
    source,umin,umax,width=electrode_indices(u,None,electrode_width_fraction)
    drain=np.where(u>=umax-width)[0]
    normal=atoms.cell[2] if layer_normal is None else layer_normal
    adj=active_adjacency(graph,normal,eta,intralayer_cutoff_A,interlayer_cutoff_A,interlayer_threshold_A)
    reached=reachable_from_targets(adj,drain)
    reachable_source=np.asarray([int(i) for i in source if int(i) in reached],dtype=int)
    return dict(source=np.asarray(source,dtype=int),drain=np.asarray(drain,dtype=int),
                reachable_source=reachable_source,reached=reached,
                source_reachable_fraction=len(reachable_source)/max(len(source),1),
                all_nodes_reachable_fraction=len(reached)/max(len(atoms),1))


def source_drain_reachability_generic(atoms, graph, direction='X',
                                      electrode_width_A=None,
                                      electrode_width_fraction=0.05):
    """Exact source-drain connectivity using every edge in a finite graph."""
    d = direction_vector(direction, atoms)
    u = np.asarray(atoms.positions) @ d
    source, umin, umax, width = electrode_indices(
        u, electrode_width_A, electrode_width_fraction
    )
    drain = np.where(u >= umax - width)[0]
    adj = [set(int(e.j) for e in edges) for edges in graph]
    reached = reachable_from_targets(adj, drain)
    reachable_source = np.asarray(
        [int(i) for i in source if int(i) in reached], dtype=int
    )
    return dict(
        source=np.asarray(source, dtype=int),
        drain=np.asarray(drain, dtype=int),
        reachable_source=reachable_source,
        reached=reached,
        source_reachable_fraction=len(reachable_source) / max(len(source), 1),
        all_nodes_reachable_fraction=len(reached) / max(len(atoms), 1),
    )
