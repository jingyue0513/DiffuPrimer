#!/usr/bin/env python3
"""Compare original multiquery and independent-query searches on actual references."""
import json
import math
import subprocess
from collections import Counter
from pathlib import Path

from blast_query_batches import fasta,run_query_batches
from diffuprimer.io import sha256,utc_now,write_json
from diffuprimer.specificity import parse_blast_tsv


root=Path(__file__).resolve().parents[1]
directory=root/'results/checks/blast_batching';directory.mkdir(parents=True,exist_ok=True)
calibration=root/'results/evaluations/primer3_calibration_validation_seed20261013_amp80-200'
queries=fasta(calibration/'specificity_queries.fasta')
selected={}
for qid,sequence in queries.items():
    counts=Counter(sequence[i:i+2] for i in range(len(sequence)-1))
    entropy=-sum(n/(len(sequence)-1)*math.log2(n/(len(sequence)-1)) for n in counts.values())
    if entropy>=3.1:selected[qid]=sequence
    if len(selected)==8:break
assert len(selected)>=4
query=directory/'queries.fasta';query.write_text(''.join(f'>{q}\n{s}\n' for q,s in selected.items()))
comparisons=[]
for name in ['transcriptome','genome']:
    whole=directory/f'{name}_multi.tsv';split=directory/f'{name}_split_{sha256(root/"scripts/blast_query_batches.py")[:8]}.tsv'
    command=[str(root/'tools/ncbi-blast-2.17.0+/bin/blastn'),'-task','blastn-short','-query',str(query),
        '-db',str(root/'data/processed/blastdb'/name),'-out',str(whole),'-dust','no','-soft_masking','false',
        '-word_size','7','-evalue','1000','-max_target_seqs','1000000','-max_hsps','1000000','-num_threads','8',
        '-outfmt','6 qseqid sseqid qstart qend sstart send length mismatch gaps qlen qseq sseq bitscore evalue']
    subprocess.run(command,check=True,capture_output=True)
    split_command=command.copy();split_command[split_command.index('-out')+1]=str(split)
    run_query_batches(split_command,validate=False,batch_size=3)
    left,right=parse_blast_tsv(whole),parse_blast_tsv(split)
    def signature(hits):
        return Counter((qid,json.dumps(row,sort_keys=True)) for qid,rows in hits.items() for row in rows)
    assert signature(left)==signature(right),(name,'Qualifying hits differ')
    comparisons.append({'reference':name,'queries':len(selected),'qualifying_hits':sum(map(len,left.values())),
        'whole_output_sha256':sha256(whole),'split_output_sha256':sha256(split),'exact_qualifying_hit_set_equality':True})
write_json(directory/'check.json',{'at':utc_now(),'exact_qualifying_hit_set_equality':True,'comparisons':comparisons,
    'script_sha256':sha256(root/'scripts/blast_query_batches.py'),'classifier_sha256':sha256(root/'src/diffuprimer/specificity.py'),
    'documentation':{'search_parameters':'https://www.ncbi.nlm.nih.gov/books/NBK279684/','threading':'https://www.ncbi.nlm.nih.gov/books/NBK571452/'}})
print('Independent queries preserve actual reference qualifying hits, scores and E values',flush=True)
