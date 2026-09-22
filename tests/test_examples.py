from pathlib import Path
import json
from hopping3d.cli import main

def test_first_passage_smoke(tmp_path):
    root=Path(__file__).resolve().parents[1]
    cfg=json.loads((root/'examples/first_passage/params.json').read_text())
    cfg['structure']['file']=str(root/'examples/first_passage/structure.cif')
    cfg['transport']['n_trajectories']=20
    cfg['transport']['voltages_V']=[0.05]
    cfg['output']['record_paths']=2
    cfg['output']['directory']=str(tmp_path/'fp')
    p=tmp_path/'fp.json'; p.write_text(json.dumps(cfg))
    assert main([str(p)])==0
    assert (tmp_path/'fp/sweep.csv').exists()

def test_diffusion_smoke(tmp_path):
    root=Path(__file__).resolve().parents[1]
    cfg=json.loads((root/'examples/diffusion_tensor/params.json').read_text())
    cfg['structure']['file']=str(root/'examples/diffusion_tensor/structure.cif')
    cfg['diffusion']['n_walkers']=20; cfg['diffusion']['observation_time_s']=1e-12
    cfg['output']['directory']=str(tmp_path/'diff')
    p=tmp_path/'diff.json'; p.write_text(json.dumps(cfg))
    assert main([str(p)])==0
    assert (tmp_path/'diff/diffusion_tensor.csv').exists()
