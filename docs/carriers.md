# Electrons and holes

Hopping3D uses a single-carrier transport model. Electron and hole transport can be requested in the same workflow, but they are propagated in **separate ensembles**.

The electronic workflow is:

```text
structure -> one EHT Hamiltonian
          -> occupied frontier projection -> hole parameters
          -> empty frontier projection    -> electron parameters
```

The transport workflow is:

```text
hole parameters     -> hole CTMC/KMC -> D_h, tau_h, paths_h, ...
electron parameters -> electron CTMC/KMC -> D_e, tau_e, paths_e, ...
```

For molecular Marcus transport, use carrier-specific `J_ij`, carrier-specific site/free-energy differences, and carrier-specific `lambda` whenever available. Ordinary EHT does not reliably determine nuclear reorganization energies.

A simultaneous interacting electron-hole simulation (recombination, excitons, carrier-carrier Coulomb terms) is outside the present model.
