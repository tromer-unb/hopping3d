# Configuration parameters

`cutoff_A`: geometric neighbor cutoff used to build the hopping graph.
`direction`: source-drain direction (`X`, `Y`, `Z`, `A`, `B`, `C`, or a vector).
`temperature_K`: temperature entering Miller-Abrahams activation.
`xi0_A`: bare localization length.
`nu0_Hz`: attempt frequency; reduced time `nu0*tau` removes its absolute scale.
`interlayer_scale`: relative transverse/interlayer hopping weight, eta_perp.
`interlayer_threshold_A`: geometric criterion distinguishing interlayer edges.
`intralayer_cutoff_A`, `interlayer_cutoff_A`: optional edge-specific cutoffs.
`site_offsets_eV`: empirical species-dependent onsite energies.
`pair_scales`: empirical relative hopping weights for chemical pairs.
`vacancy_fraction`: fraction of atoms removed before graph construction.
`n_trajectories`: number of finite-device KMC trajectories.
`max_time_s`, `max_steps`: trajectory stopping criteria.

The model is phenomenological: `interlayer_scale`, site offsets, and pair scales are model parameters, not claimed DFT transfer integrals.
