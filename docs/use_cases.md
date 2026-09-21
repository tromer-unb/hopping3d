# Choosing the physical model

This page answers the question that should be asked before running any KMC calculation:

> **What is a hopping site in my system, and what physical quantity should determine the rate between two sites?**

The numerical engine can propagate trajectories on many graphs, but a graph is only physically meaningful if its nodes and edge rates correspond to the localized states relevant to the material.

## Decision table

| Question | Recommended calculation | Site representation | Notes |
|---|---|---|---|
| Periodic 3D bulk transport | `diffusion_tensor` | usually atomic or pre-coarse-grained sites | returns the full diffusion tensor and principal axes |
| Finite source-to-drain device | `first_passage` | atomic or pre-coarse-grained sites | returns passage fraction, first-passage time and path statistics |
| Is the graph connected at all? | `reachability` | same graph used by the transport model | no stochastic KMC time scale is required |
| Layered atomic material | any of the above | atoms | use `interlayer_scale` and layer geometry when appropriate |
| Atomic impurity in a host, e.g. Li on graphene | usually `first_passage` or `diffusion_tensor` | host atoms plus an explicit impurity model | decide whether the impurity is a transport site, an energy perturbation, or both |
| Molecular crystal / amorphous molecular box | typically `diffusion_tensor` | usually **one molecule = one site** | molecular coarse graining and orientation-dependent couplings are not automatic in v1.0 |
| Finite curved object, e.g. nanocone | `first_passage` or `reachability` | atoms | no planar geometry is required |

## 1. Atomic networks

The default public model uses

```text
one atom = one hopping site
```

Neighbor edges are defined from geometry, and the default rate is Miller--Abrahams-like. This is the representation used for the atomic-network calculations associated with the repository.

This mode is directly usable for graphene-like networks, layered atomic crystals, defected atomic lattices, finite nanostructures, and other situations where localized states are intentionally associated with atomic positions.

## 2. Molecular systems: a benzene box

For charge transfer through a molecular solid, the most natural coarse graining is commonly

```text
one molecule = one hopping site
```

rather than one carbon atom = one hopping site.

A benzene box containing 500 molecules would therefore contain approximately 500 transport sites after molecular coarse graining, not 6000 atomic sites.

### Why molecular orientation matters

Two benzene molecules with the same center-to-center distance should generally **not** be assumed to have the same electronic coupling if one pair is face-to-face and another is edge-on or nearly perpendicular.

A purely radial model,

```text
J_ij = J(r_ij)
```

cannot distinguish these configurations.

A better molecular model uses

```text
J_ij = J(r_ij, orientation_i, orientation_j, displacement_direction)
```

or supplies a transfer integral `J_ij` calculated externally for each molecular pair.

For planar aromatic molecules, useful geometric descriptors include

- center-to-center distance `r_ij`;
- angle between molecular plane normals `n_i . n_j`;
- slip / lateral displacement between ring centers;
- angle between the intermolecular displacement and each molecular normal.

A simple empirical orientation factor could be used for exploratory work, but quantitative organic-semiconductor calculations should preferably use electronic couplings from a suitable electronic-structure or fragment-orbital model.

### Miller--Abrahams versus Marcus

The current public engine uses a Miller--Abrahams-type rate. Molecular charge transfer is often represented instead with a Marcus-type rate,

```text
k_ij = (2*pi/hbar) |J_ij|^2 / sqrt(4*pi*lambda*kB*T)
       * exp[-(DeltaG_ij + lambda)^2 / (4*lambda*kB*T)]
```

which introduces

- the pair transfer integral `J_ij`;
- the molecular reorganization energy `lambda`;
- the site-energy difference / driving force `DeltaG_ij`.

The present repository does not claim that a raw atomistic benzene CIF plus the default atomic rate is a quantitative benzene transport model. A molecular coarse-grained layer with orientation-dependent or externally supplied `J_ij` is the correct extension.

## 3. Graphene with Li impurities or adsorbates

Several distinct physical models are possible, and they should not be mixed silently.

### Model A: Li is only an electrostatic / chemical perturbation

Carbon atoms remain the hopping sites. The Li atom modifies nearby carbon site energies,

```text
Li position -> local potential on neighboring C sites -> modified C-to-C rates
```

This is often the cleanest reduced model when Li acts primarily as a dopant or local perturbation rather than as a transport intermediate.

The current engine can represent this with site-energy corrections, including per-species, per-index, or explicitly supplied arrays.

### Model B: Li is itself a hopping site

Both C and Li are nodes,

```text
C <-> C
C <-> Li
Li <-> Li
```

and the pair-dependent relative couplings can differ through `pair_scales`.

For example, one may use independent phenomenological weights for `C-C`, `C-Li`, and `Li-Li`. Such weights must be physically justified or treated explicitly as sensitivity parameters.

### Model C: Li is a perturbation and a transport intermediate

This is also possible: Li can have its own site energy while simultaneously shifting nearby C site energies. This is the most flexible reduced model but introduces more parameters.

For a quantitative prediction, site energies and C--Li couplings should ideally be calibrated by an electronic-structure calculation or experiment.

## 4. Finite curved systems: nanocones, cages, nanoparticles

The finite-device calculation does not require a flat sheet or a periodic cell.

A carbon nanocone can be represented directly as a finite atomic graph:

```text
nanocone geometry -> finite neighbor graph -> source/drain projection -> KMC first passage
```

Choose a transport direction, for example the cone axis. Source and drain are then defined from the extrema of

```text
u_i = r_i . n_transport
```

where `n_transport` is the selected direction vector.

This makes it possible to study, for example,

- transport from the cone apex to the base;
- transport across the base diameter;
- connectivity loss after removing atoms;
- path localization around topological defects;
- direction-dependent first-passage statistics.

For a finite object use `first_passage` or `reachability`, not `diffusion_tensor`. The periodic diffusion tensor is intended for a periodically repeated bulk system.

## 5. What the program does not decide for you

The program can construct and propagate a graph, but it cannot infer the correct electronic coarse graining from geometry alone.

The user must decide:

1. what constitutes a localized transport site;
2. which pairs are allowed to exchange charge;
3. how site energies are assigned;
4. what determines the pair rate;
5. whether the system is finite or periodic;
6. whether the desired observable is connectivity, first passage, or bulk diffusion.

These choices are part of the physical model, not merely software settings.
