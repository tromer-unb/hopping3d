# Configuration parameters

A normal user runs

```bash
hopping3d params.json
```

The JSON file selects one of three public calculation modes.

## Top-level calculation modes

`calculation: "first_passage"`
: finite source-to-drain KMC. Use for devices, ribbons, flakes, finite clusters, and finite 3D networks.

`calculation: "diffusion_tensor"`
: periodic bulk KMC. Returns the full diffusion tensor and principal diffusivities/axes.

`calculation: "reachability"`
: graph connectivity only. No KMC time scale is sampled.

## Structure

`structure.file`
: preferred path to the structure. ASE-readable formats such as CIF, XYZ and POSCAR are supported when appropriate for the selected calculation.

`structure.cif`
: legacy alias retained for backward compatibility.

`structure.repeat`
: cell replication factors. Primarily useful for periodic/crystalline structures.

`structure.vacancy_fraction`
: fraction of candidate atoms removed before graph construction.

`structure.vacancy_seed`
: random seed for reproducible vacancy placement.

## Graph

`cutoff_A`
: geometric neighbor cutoff used to build the hopping graph.

`interlayer_scale`
: relative transverse/interlayer hopping weight, `eta_inter`.

`layer_normal`
: vector used to classify transverse edges in layered models.

`interlayer_threshold_A`
: geometric criterion distinguishing interlayer edges.

`intralayer_cutoff_A`, `interlayer_cutoff_A`
: optional edge-specific cutoffs.

## Transport

`direction`
: source-drain direction (`X`, `Y`, `Z`, `A`, `B`, `C`, or a 3-vector).

`temperature_K`
: temperature entering Miller--Abrahams activation.

`xi0_A`
: bare localization length.

`nu0_Hz`
: attempt frequency; reduced time `nu0*tau` removes its global clock scale.

`pair_scales`
: phenomenological relative hopping weights for chemical pairs.

`n_trajectories`
: number of finite-device KMC trajectories.

`max_time_s`, `max_steps`
: finite-trajectory stopping criteria.

`electrode_width_A` or `electrode_width_fraction`
: source/drain slab width for finite-device calculations.

## Chemistry

`site_offsets_eV`
: empirical species-dependent onsite energies.

`site_offsets_by_index_eV`
: energy corrections for specific atomic-site indices.

`site_offsets_array_eV`
: full site-energy correction array supplied by the user.

### Explicit EHT / external couplings

`edge_transfer_integrals_eV`
: mapping from an edge key such as `"12-27"` to an explicit transfer integral. When supplied, the rate weight is proportional to `|J_ij/J_ref|^2`.

`J_ref_eV`
: reference coupling used to make pair weights dimensionless for Miller--Abrahams calculations.

`strict_edge_transfer_integrals`
: if true, graph edges absent from the explicit coupling map receive zero electronic coupling rather than a default unit weight.

`rate_model`
: `"miller_abrahams"` or `"marcus"`.

`reorganization_energy_eV`
: Marcus reorganization energy. This is not predicted reliably by ordinary EHT and should be supplied externally.

`edge_reorganization_energy_eV`
: optional edge-specific Marcus reorganization energies.

The EHT layer in `hopping3d_eht` can infer atomic and molecular couplings for the validated article use cases. See `docs/eht_to_transport.md` and `docs/use_cases.md`. Raw EHT level differences are screening quantities unless independently validated.
