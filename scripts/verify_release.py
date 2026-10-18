#!/usr/bin/env python3
"""Offline verification of release file hashes and candidate/screen coverage."""
import hashlib,json,zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def digest(data):return hashlib.sha256(data).hexdigest()

def main():
    manifest=json.loads((ROOT/'release/data_manifest.json').read_text())
    checked=0;screens=0;candidates=0
    for record in manifest['archives']:
        path=ROOT/'release'/record['archive'];assert digest(path.read_bytes())==record['sha256'],path
        with zipfile.ZipFile(path) as z:
            expected={x['path']:x for x in record['members']};assert set(z.namelist())==set(expected)
            for name,item in expected.items():
                data=z.read(name)
                assert len(data)==item['bytes'] and digest(data)==item['sha256'],name
                checked+=1
                if name.endswith('/specificity_bowtie_status.json'):
                    status=json.loads(data);base=name.rsplit('/',1)[0]
                    raw=z.read(base+'/candidates.jsonl');scored=z.read(base+'/specificity_bowtie.jsonl')
                    assert status['status']=='completed',name
                    assert digest(raw)==status['candidate_sha256'],name
                    assert digest(scored)==status['output_sha256'],name
                    a=[json.loads(x) for x in raw.splitlines()];b=[json.loads(x) for x in scored.splitlines()]
                    assert {x['candidate_id'] for x in a}=={x['candidate_id'] for x in b},name
                    screens+=1;candidates+=len(a)
    result={'passed':True,'archives':len(manifest['archives']),'files':checked,
            'completed_screen_runs_verified':screens,'screened_candidate_records':candidates}
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
