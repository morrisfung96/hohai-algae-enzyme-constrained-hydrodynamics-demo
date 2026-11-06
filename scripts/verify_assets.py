"""Verify the published file inventory and excluded-model boundary."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda: f.read(1 << 20), b''):
            h.update(part)
    return h.hexdigest()

def main():
    manifest_path = ROOT / 'manifests/files.json'
    if not manifest_path.is_file():
        raise RuntimeError('Release file inventory is missing; this directory is not a finalized package.')
    records = json.loads(manifest_path.read_text(encoding='utf-8'))['files']
    errors = []
    for r in records:
        p = ROOT / r['path']
        if not p.is_file():
            errors.append({'path': r['path'], 'error': 'missing'})
        elif p.stat().st_size != r['bytes'] or sha256(p) != r['sha256']:
            errors.append({'path': r['path'], 'error': 'size or checksum differs'})
    excluded = [str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file()
                and p.suffix.lower() in {'.xml', '.sbml', '.ipynb', '.pkl', '.pickle'}
                and not any(x in p.parts for x in ['.git', '.venv', 'external_data'])]
    report = {'status': 'PASS' if not errors and not excluded else 'FAIL',
              'inventory_files': len(records), 'errors': errors,
              'excluded_format_files': excluded}
    print(json.dumps(report, indent=2))
    if errors or excluded:
        raise SystemExit(1)

if __name__ == '__main__':
    main()
