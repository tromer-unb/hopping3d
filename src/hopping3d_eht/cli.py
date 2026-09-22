"""Direct Extended-Huckel parametrization CLI for Hopping3D."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from ase.io import read
from ase.build import molecule
from .alignment import molecule_graphene_alignment_eht
from .interface import parameterize_interface,render_interface_markdown
from .trajectory import process_trajectory,save_trajectory_report
from .multi_interface import parameterize_coverage_interface,render_coverage_markdown

def _dump(obj,output=None):
    text=json.dumps(obj,indent=2)
    if output:
        Path(output).write_text(text+'\n'); print(f'Wrote {output}')
    else: print(text)

def cmd_doctor(_):
    def probe(names):
        out={}
        for name in names:
            try:
                m=__import__(name); out[name]=getattr(m,'__version__','installed')
            except Exception as exc: out[name]=f'MISSING: {exc}'
        return out
    required=probe(('numpy','scipy','ase')); optional=probe(('pandas','matplotlib')); ready=all(not str(v).startswith('MISSING') for v in required.values())
    _dump({'python':sys.version.split()[0],'required_python_dependencies':required,'optional_python_dependencies':optional,'portable_interface_ready':ready,
           'electronic_backend':'direct Extended-Huckel','external_executables_required':[],'basic_runtime':'EHT-only interface/coverage/trajectory require NumPy + SciPy + ASE'})

def cmd_graphene(ns):
    from .parameterize import parameterize_graphene
    out=parameterize_graphene(read(ns.structure)); out['user_mode']='portable_eht'; out['external_executables_required']=[]; _dump(out,ns.output)

def cmd_align(ns):
    g=read(ns.graphene); m=read(ns.molecule) if ns.molecule else molecule('C6H6'); m.positions-=m.positions.mean(axis=0); m.pbc=False
    out=molecule_graphene_alignment_eht(g,m); out['user_mode']='portable_eht_with_shipped_offline_calibration'; out['external_executables_required']=[]; _dump(out,ns.output)

def cmd_interface(ns):
    out=parameterize_interface(read(ns.structure),patch_radius_A=ns.patch_radius); out['input_structure']=str(ns.structure); _dump(out,ns.output)
    if ns.report: Path(ns.report).write_text(render_interface_markdown(out)); print(f'Wrote {ns.report}')

def cmd_coverage(ns):
    out=parameterize_coverage_interface(read(ns.structure),patch_radius_A=ns.patch_radius,electron_offset_eV=ns.electron_offset); out['input_structure']=str(ns.structure); _dump(out,ns.output)
    if ns.report: Path(ns.report).write_text(render_coverage_markdown(out)); print(f'Wrote {ns.report}')

def cmd_trajectory(ns):
    out=process_trajectory(ns.trajectory,stride=ns.stride,max_frames=ns.max_frames,patch_radius_A=ns.patch_radius,electron_offset_eV=ns.electron_offset); save_trajectory_report(out,ns.output,ns.csv); print(f'Wrote {ns.output}')
    if ns.csv: print(f'Wrote {ns.csv}')

def main(argv=None):
    ap=argparse.ArgumentParser(prog='hopping3d-eht',description='Direct Extended-Huckel parametrization for Hopping3D transport'); sub=ap.add_subparsers(dest='command',required=True)
    p=sub.add_parser('doctor',help='check the EHT runtime'); p.set_defaults(func=cmd_doctor)
    p=sub.add_parser('graphene',help='estimate graphene EHT hopping parameters from a structure'); p.add_argument('structure'); p.add_argument('-o','--output',default='hopping3d_graphene_eht.json'); p.set_defaults(func=cmd_graphene)
    p=sub.add_parser('align',help='rough conjugated-molecule/graphene level alignment'); p.add_argument('graphene'); p.add_argument('--molecule'); p.add_argument('-o','--output',default='hopping3d_alignment.json'); p.set_defaults(func=cmd_align)
    p=sub.add_parser('interface',help='auto-detect and parameterize a graphene + molecular adsorbate interface'); p.add_argument('structure'); p.add_argument('--patch-radius',type=float,default=8.0); p.add_argument('-o','--output',default='hopping3d_interface.json'); p.add_argument('--report'); p.set_defaults(func=cmd_interface)
    p=sub.add_parser('coverage',help='parameterize graphene with multiple molecular adsorbates'); p.add_argument('structure'); p.add_argument('--patch-radius',type=float,default=8.0); p.add_argument('--electron-offset',type=float); p.add_argument('-o','--output',default='hopping3d_coverage.json'); p.add_argument('--report'); p.set_defaults(func=cmd_coverage)
    p=sub.add_parser('trajectory',help='process a multi-frame molecular-interface trajectory with direct EHT'); p.add_argument('trajectory'); p.add_argument('--stride',type=int,default=1); p.add_argument('--max-frames',type=int); p.add_argument('--patch-radius',type=float,default=8.0); p.add_argument('--electron-offset',type=float); p.add_argument('-o','--output',default='hopping3d_trajectory.json'); p.add_argument('--csv'); p.set_defaults(func=cmd_trajectory)
    ns=ap.parse_args(argv)
    try: ns.func(ns)
    except ValueError as exc: ap.error(str(exc))

if __name__=='__main__': main()
