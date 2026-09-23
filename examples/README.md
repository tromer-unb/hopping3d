# Examples

The repository examples correspond directly to the systems discussed in the transport framework.

- `first_passage/`: finite source-to-drain KMC.
- `diffusion_tensor/`: periodic diffusion tensor.
- `reachability/`: graph connectivity without a kinetic time scale.
- `article_systems/carbon_stack/`: layered carbon and dimensional rescue.
- `article_systems/bn/`: 3D BN tensor benchmark.
- `article_systems/w2o6/`: chemically heterogeneous W-O network.
- `article_systems/benzene_box/`: generated 100-benzene periodic molecular benchmark.
- `article_systems/graphene_benzene/`: one or several benzene molecules on graphene.

The 100-benzene periodic box is generated programmatically by `reproduction/article_cases/run_benzene_box.py`, so no large pre-generated molecular box is required in the repository.

The most complete reproducibility entry point is `reproduction/article_cases/run_all.sh`.
