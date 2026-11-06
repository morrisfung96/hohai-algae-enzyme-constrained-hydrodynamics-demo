"""Verify and unpack a published yearly input archive; no network required."""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda:f.read(1 << 20),b''): h.update(part)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--output', type=Path, default=ROOT / 'external_data')
    args = parser.parse_args()
    inventory = json.loads((ROOT / 'manifests/external_inputs.json').read_text())
    records = [r for r in inventory['assets'] if r['file'] == args.archive.name]
    if len(records) != 1: parser.error('Archive name is not part of the formal input inventory.')
    expected = records[0]
    if args.archive.stat().st_size != expected['bytes'] or sha256(args.archive) != expected['sha256']:
        raise RuntimeError('Archive checksum or size differs from the release manifest.')
    dest = args.output.resolve()
    with ZipFile(args.archive) as archive:
        for item in archive.infolist():
            target = (dest / item.filename).resolve()
            if not target.is_relative_to(dest): raise RuntimeError('Unsafe path in input archive.')
        archive.extractall(dest)
    folder = dest / str(expected['year'])
    for record in json.loads((folder / 'files.json').read_text())['files']:
        p = folder / record['path']
        if not p.is_file() or p.stat().st_size != record['bytes'] or sha256(p) != record['sha256']:
            raise RuntimeError('Extracted input differs: ' + record['path'])
    print(json.dumps({'status':'PASS','year':expected['year'],'verified_extracted_inputs':True}))

if __name__ == '__main__':
    main()
