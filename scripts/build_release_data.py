#!/usr/bin/env python3
"""Package frozen data/evidence without downloaded references, caches or vendors."""
import hashlib
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'release'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def archive(name, paths):
    paths=sorted({p for p in paths if p.is_file() and '__pycache__' not in p.parts})
    rows=[]
    with zipfile.ZipFile(OUT/name,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in paths:
            rel=str(p.relative_to(ROOT));z.write(p,rel);rows.append({'path':rel,'bytes':p.stat().st_size,'sha256':sha(p)})
    return {'archive':name,'bytes':(OUT/name).stat().st_size,'sha256':sha(OUT/name),'members':rows}


def main():
    OUT.mkdir(exist_ok=True)
    records=[]
    # Data restore to their historical paths so frozen hashes stay meaningful.
    names=['targets_pilot.jsonl','attention_external_targets.jsonl','dataset_partial_pilot.jsonl',
           'dataset_partial_pilot.manifest.json','pilot_gene_selection.jsonl','pilot_gene_selection.manifest.json',
           'reference_manifest.json','acquisition_ledger.jsonl','primerbank_pilot.jsonl','pilot_all_isoforms.fasta',
           'computational_candidates.manifest.json','computational_candidates_ledger.jsonl']
    data=[ROOT/'data/processed'/n for n in names]
    data+=list((ROOT/'data/splits').glob('*'))
    data+=[ROOT/f'results/checkpoints/{m}_weak_partial_pilot_seed20261013.pt' for m in ['direct','diffusion']]
    records.append(archive('benchmark_data.zip',data))
    for study,zipname in [('attention','original_evidence.zip'),('attention_external','additional_evidence.zip')]:
        base=ROOT/'experiments'/study
        summary=base/'results/study_summary.json';s=json.loads(summary.read_text())
        paths=[summary]
        paths += list((base/'results').glob('study_summary.*'))
        # Every candidate source used in the summaries, not tool-native working files.
        dirs={ROOT/str(g['run']).split(str(ROOT)+'/')[-1] for g in s['summaries']}
        for k in ['inputs','validation_inputs','evaluation_directories']:
            value=s.get(k,[])
            if isinstance(value,list):
                for x in value:
                    if isinstance(x,str) and 'evaluations/' in x:
                        dirs.add(ROOT/x.split(str(ROOT)+'/')[-1])
        for d in dirs:
            for name in ['manifest.json','completion.json','candidates.jsonl','per_gene.jsonl',
                         'specificity_bowtie.jsonl','specificity_bowtie_status.json']:
                paths.append(d/name)
        paths+=list((base/'results/exon_interpretation').glob('*.json'))
        paths+=list((base/'results/exon_interpretation').glob('*.jsonl'))
        records.append(archive(zipname,paths))
    mech=ROOT/'experiments/mechanism'
    paths=[p for p in mech.rglob('*') if p.suffix in {'.json','.jsonl','.csv','.txt','.pdf','.png'}]
    records.append(archive('mechanism_seed42.zip',paths))
    (OUT/'data_manifest.json').write_text(json.dumps({'format_version':1,'archives':records},indent=2)+'\n')
    for r in records:print(r['archive'],r['bytes'],len(r['members']),flush=True)


if __name__=='__main__':main()
