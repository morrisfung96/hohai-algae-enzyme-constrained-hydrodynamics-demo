"""Assemble the approved artwork; this is not a new numerical renderer."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/figures')
    args = parser.parse_args()
    records = {r['path']: r for r in json.loads((ROOT / 'manifests/asset_sources.json').read_text())}
    args.output.mkdir(parents=True, exist_ok=True)
    out = []
    for name in ['Figure1.png', 'Figure2.png', 'Figure3.png', 'Figure4.png', 'FigureS6_1.png', 'FigureS6_2.png']:
        relative = 'artwork/approved/' + name
        src = ROOT / relative
        sha = hashlib.sha256(src.read_bytes()).hexdigest()
        if sha != records[relative]['sha256']:
            raise RuntimeError('Artwork checksum mismatch: ' + relative)
        shutil.copyfile(src, args.output / name)
        out.append({'file': name, 'sha256': sha, 'operation': 'approved artwork assembly'})
    (args.output / 'figure_reproduction.json').write_text(json.dumps({
        'status': 'ASSEMBLED_FROM_APPROVED_ARTWORK',
        'numerical_model_or_LP_run': False, 'files': out,
        'note': 'Approved schematic/network artwork is retained. This command does not claim to recompute every quantitative panel from model states.'
    }, indent=2) + '\n')
    print(json.dumps({'status': 'ASSEMBLED_FROM_APPROVED_ARTWORK', 'figure_files': len(out)}))

if __name__ == '__main__':
    main()
