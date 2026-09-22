"""EHT-only processing of multi-frame molecular-interface trajectories."""
from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
from ase.io import read
from .multi_interface import parameterize_coverage_interface

def _stats(values):
    x=np.asarray([v for v in values if v is not None],float)
    if len(x)==0: return {'n':0,'mean':None,'std':None,'min':None,'max':None}
    return {'n':int(len(x)),'mean':float(x.mean()),'std':float(x.std()),'min':float(x.min()),'max':float(x.max())}

def summarize_frame(result,frame_index):
    gm=[a.get('parameters',{}).get('electron_J_eff_eV') for a in result['adsorbates']]; mm=[p.get('electron_J_eff_eV') for p in result['molecule_pairs']]
    row={'frame':int(frame_index),'n_adsorbates':int(result['detection']['n_adsorbates']),'coverage_nm-2':result.get('coverage',{}).get('molecules_per_nm2'),'GM_J_eV':_stats(gm),'MM_J_eV':_stats(mm)}
    dyn=result.get('electron_coverage_dynamics')
    if dyn:
        m=dyn['metrics']; row['electron_offset_eV']=dyn['energy_offset_eV']; row['dynamics']={k:m[k] for k in ('total_capture_rate_s-1','mean_wait_to_capture_s','mean_molecular_excursion_s','probability_any_inter_molecular_hop','expected_inter_molecular_hops_per_excursion','probability_return_via_different_molecule','rms_lateral_migration_A','molecular_stationary_probability')}
    return row

def process_trajectory(path,stride=1,max_frames=None,patch_radius_A=8.0,electron_offset_eV=None):
    frames=read(str(path),index=':'); frames=frames if isinstance(frames,list) else [frames]; sel=frames[::max(int(stride),1)]
    if max_frames is not None: sel=sel[:int(max_frames)]
    out=[]
    for k,atoms in enumerate(sel):
        r=parameterize_coverage_interface(atoms,patch_radius_A=patch_radius_A,electron_offset_eV=electron_offset_eV); out.append(summarize_frame(r,k))
    return {'schema':'hopping3d-trajectory-0.1','input':str(path),'n_input_frames':int(len(frames)),'n_processed_frames':int(len(out)),'stride':int(stride),'electron_offset_eV':electron_offset_eV,'electronic_parameterization':'direct EHT','frames':out}

def aggregate_trajectory(report):
    rows=report['frames']; gm=np.asarray([r['GM_J_eV']['mean'] for r in rows if r['GM_J_eV']['mean'] is not None],float); mm=np.asarray([r['MM_J_eV']['mean'] for r in rows if r['MM_J_eV']['mean'] is not None],float)
    out={'n_frames':len(rows),'GM_J_mean_eV':float(gm.mean()) if len(gm) else None,'GM_J_std_over_frames_eV':float(gm.std()) if len(gm) else None,'MM_J_mean_eV':float(mm.mean()) if len(mm) else None,'MM_J_std_over_frames_eV':float(mm.std()) if len(mm) else None}
    if rows and 'dynamics' in rows[0]:
        keys=rows[0]['dynamics'].keys(); out['dynamics_mean']={k:float(np.mean([r['dynamics'][k] for r in rows])) for k in keys}; out['dynamics_std']={k:float(np.std([r['dynamics'][k] for r in rows])) for k in keys}
    return out

def write_trajectory_csv(report,path):
    flat=[]
    for r in report['frames']:
        q={'frame':r['frame'],'n_adsorbates':r['n_adsorbates'],'coverage_nm-2':r['coverage_nm-2'],'GM_J_mean_eV':r['GM_J_eV']['mean'],'GM_J_std_eV':r['GM_J_eV']['std'],'MM_J_mean_eV':r['MM_J_eV']['mean'],'MM_J_std_eV':r['MM_J_eV']['std']}
        if 'dynamics' in r: q.update(r['dynamics'])
        flat.append(q)
    fields=list(flat[0].keys()) if flat else ['frame']
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(flat)

def save_trajectory_report(report,json_path,csv_path=None):
    report=dict(report); report['aggregate']=aggregate_trajectory(report); Path(json_path).write_text(json.dumps(report,indent=2)+'\n')
    if csv_path: write_trajectory_csv(report,csv_path)
    return report
