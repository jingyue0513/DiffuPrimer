#!/usr/bin/env python3
"""Verify and restore packaged benchmark data; refuse changed existing files."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]

def sha(data):return hashlib.sha256(data).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT,help='Destination checkout (default: repository root)')
    a=p.parse_args();root=a.root.resolve()
    manifest=json.loads((ROOT/'release/data_manifest.json').read_text())
    count=0
    for record in manifest['archives']:
        archive=ROOT/'release'/record['archive']
        assert sha(archive.read_bytes())==record['sha256'],archive
        with zipfile.ZipFile(archive) as z:
            expected={r['path']:r for r in record['members']}
            assert set(z.namelist())==set(expected)
            for name,item in expected.items():
                target=(root/name).resolve()
                if not target.is_relative_to(root):raise ValueError('Unsafe archive path')
                data=z.read(name)
                assert len(data)==item['bytes'] and sha(data)==item['sha256'],name
                if target.exists():
                    if sha(target.read_bytes())!=item['sha256']:raise ValueError(f'Existing file differs: {target}')
                else:
                    target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
                count+=1
    print(f'Verified/restored {count} frozen files in {root}')

if __name__=='__main__':main()
