#!/usr/bin/env python3
"""Freeze new RefSeq genes independently of all primer-design outcomes."""
import argparse
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT / "src"))
from diffuprimer.acquisition import load_reference_sequences, reference_gene_catalog
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl


def write_fasta(path, accessions, sequences):
    path.write_text("".join(f">{a}\n{sequences[a]}\n" for a in sorted(accessions)))


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--seed",type=int,default=20261013)
    p.add_argument("--genes",type=int,default=75)
    p.add_argument("--batch-size",type=int,default=128)
    p.add_argument("--threads",type=int,default=8)
    args=p.parse_args()
    out=Path(__file__).resolve().parent
    destination=ROOT / "data/processed/attention_external_targets.jsonl"
    source_dir=ROOT / "data/raw/refseq"
    gff=source_dir / "GCF_000001635.27_GRCm39_genomic.gff.gz"
    rna=source_dir / "GCF_000001635.27_GRCm39_rna.fna.gz"
    original_path=ROOT / "data/processed/pilot_gene_selection.jsonl"
    split_path=ROOT / "data/splits/pilot.jsonl"
    original=read_jsonl(original_path)
    original_genes={g["gene_id"] for g in original}
    train_genes={g["gene_id"] for g in read_jsonl(split_path) if g["split"]=="train"}
    train_accessions={a:g["gene_id"] for g in original if g["gene_id"] in train_genes for a in g["reference_transcripts"]}
    reference_verified=json.loads((ROOT / "results/reference_verification.json").read_text())
    expected={f["name"]:f for f in reference_verified["files"]}
    reference_hashes={f.name:sha256(f) for f in [gff,rna]}
    assert all(reference_hashes[name]==expected[name]["sha256"] for name in reference_hashes)
    catalog=reference_gene_catalog(gff)
    candidates=[g for g in catalog if g["gene_id"] not in original_genes]
    random.Random(args.seed).shuffle(candidates)
    protocol={"created_at":utc_now(),"seed":args.seed,"requested_genes":args.genes,
        "selection_method":"Shuffle the full gene-ID-sorted official protein-coding NM_ catalog after excluding all original 500 gene IDs, then take the first eligible 75 genes in that random order",
        "eligibility":"Complete official NM_ isoform sequence availability and no MMseqs2 >=80% identity, >=80% query coverage AND >=80% target coverage hit to an original training-gene isoform",
        "template_rule":"Lexicographically first annotated NM_ accession version per gene, as in the original study",
        "candidate_outcomes_consulted":False,"model_evaluation_performed":False,
        "original_selection_sha256":sha256(original_path),"original_split_sha256":sha256(split_path),
        "reference_files":[expected[gff.name],expected[rna.name]],"current_reference_sha256":reference_hashes,
        "catalog_genes":len(catalog),"candidate_genes_after_original_ID_exclusion":len(candidates),
        "original_genes":len(original_genes),"original_training_genes":len(train_genes),
        "homology_rule":{"identity":0.8,"query_coverage":0.8,"target_coverage":0.8,"isoforms":"all annotated NM_ isoforms","search_directions":"new-to-train and train-to-new"}}
    write_json(out / "protocol.json",protocol)
    write_jsonl(out / "candidate_order.jsonl",[{**g,"random_order_index":i} for i,g in enumerate(candidates)])
    train_sequences=load_reference_sequences(rna,set(train_accessions))
    assert set(train_sequences)==set(train_accessions)
    train_fasta=out / "original_training_all_isoforms.fasta"
    write_fasta(train_fasta,train_accessions,train_sequences)
    tool=ROOT / "tools/mmseqs/bin/mmseqs"
    selected=[]; ledger=[]; rejected=[]; all_edges=[]; commands=[]; all_selected_sequences={}
    for offset in range(0,len(candidates),args.batch_size):
        batch=candidates[offset:offset+args.batch_size]
        accession_gene={a:g["gene_id"] for g in batch for a in g["reference_transcripts"]}
        sequences=load_reference_sequences(rna,set(accession_gene))
        available=[g for g in batch if all(a in sequences for a in g["reference_transcripts"])]
        available_accessions={a for g in available for a in g["reference_transcripts"]}
        batch_dir=out / "homology" / f"batch_{offset:05d}"
        batch_dir.mkdir(parents=True,exist_ok=True)
        query=batch_dir / "new_all_isoforms.fasta"
        write_fasta(query,available_accessions,sequences)
        excluded=set(); batch_edges=[]
        for direction,query_fasta,target_fasta in [("new_to_train",query,train_fasta),("train_to_new",train_fasta,query)]:
            hits=batch_dir / f"{direction}.tsv"
            command=[str(tool),"easy-search",str(query_fasta),str(target_fasta),str(hits),str(batch_dir / f"tmp_{direction}"),
                "--search-type","3","--min-seq-id","0.8","-c","0.8","--cov-mode","0",
                "--max-seqs","100000","--threads",str(args.threads),"--format-output","query,target,fident,alnlen,qcov,tcov"]
            with (batch_dir / f"{direction}.log").open("w") as log:
                subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
            commands.append({"command":command,"hits_sha256":sha256(hits)})
            for line in hits.read_text().splitlines():
                q,t,identity,length,qcov,tcov=line.split("\t")
                if min(float(identity),float(qcov),float(tcov)) < 0.8:
                    continue
                new_acc,train_acc=(q,t) if direction=="new_to_train" else (t,q)
                gene=accession_gene[new_acc];excluded.add(gene)
                batch_edges.append({"new_gene":gene,"training_gene":train_accessions[train_acc],
                    "new_isoform":new_acc,"training_isoform":train_acc,"direction":direction,
                    "identity":float(identity),"query_coverage":float(qcov),"target_coverage":float(tcov),"alignment_length":int(length)})
        all_edges.extend(batch_edges)
        for local,gene in enumerate(batch):
            gid=gene["gene_id"];available_gene=all(a in sequences for a in gene["reference_transcripts"])
            reason="missing_annotated_RNA" if not available_gene else ("training_homology" if gid in excluded else "eligible")
            entry={"gene_id":gid,"random_order_index":offset+local,"source_sequence_complete":available_gene,
                "training_homology_detected":gid in excluded,"status":reason}
            if reason=="eligible" and len(selected)<args.genes:
                accession=gene["reference_transcripts"][0]
                record={**gene,"transcript_accession_version":accession,"template_sequence":sequences[accession],
                    "template_sha256":hashlib.sha256(sequences[accession].encode()).hexdigest(),
                    "split":"external_test","cluster_id":"external_"+gid,"external_random_order_index":offset+local,
                    "target_selection_rule":protocol["template_rule"],"reference_assembly":"GCF_000001635.27",
                    "annotation_version":"GCF_000001635.27-RS_2024_02"}
                selected.append(record)
                all_selected_sequences.update({a:sequences[a] for a in gene["reference_transcripts"]})
                entry["selected"]=True
            else:
                entry["selected"]=False
                if reason!="eligible":rejected.append(entry)
            ledger.append(entry)
        write_jsonl(out / "selection_ledger.jsonl",ledger)
        print(json.dumps({"examined_genes":len(ledger),"selected_genes":len(selected),"rejected_homology":sum(r["status"]=="training_homology" for r in ledger)}),flush=True)
        if len(selected)==args.genes:break
    if len(selected)!=args.genes:
        raise ValueError("Requested cohort size unavailable in the complete source catalog")
    write_jsonl(destination,selected)
    write_jsonl(out / "gene_selection.jsonl",[{k:v for k,v in row.items() if k!="template_sequence"} for row in selected])
    write_jsonl(out / "external_split.jsonl",[{"gene_id":r["gene_id"],"split":"external_test","cluster_id":r["cluster_id"]} for r in selected])
    write_jsonl(out / "homology_edges.jsonl",all_edges)
    selected_fasta=out / "selected_all_isoforms.fasta"
    write_fasta(selected_fasta,all_selected_sequences,all_selected_sequences)
    selected_genes={r["gene_id"] for r in selected}
    checks={"targets":len(selected),"gene_ID_overlap_with_original_500":len(selected_genes&original_genes),
        "gene_ID_overlap_with_original_training":len(selected_genes&train_genes),
        "selected_training_homology_edges":sum(r["new_gene"] in selected_genes for r in all_edges),
        "all_selected_isoforms_available":all(a in all_selected_sequences for r in selected for a in r["reference_transcripts"]),
        "all_templates_follow_lexicographic_rule":all(r["transcript_accession_version"]==min(r["reference_transcripts"]) for r in selected),
        "all_reference_hashes_match_official_verified_snapshot":True,
        "primer_feasibility_or_model_outputs_used_in_selection":False,"model_evaluation_performed":False}
    assert checks["gene_ID_overlap_with_original_500"]==0 and checks["selected_training_homology_edges"]==0
    write_json(out / "source_check.json",checks)
    manifest={**protocol,"completed_at":utc_now(),"status":"frozen_external_gene_cohort",
        "targets_path":str(destination),"targets_sha256":sha256(destination),"selected_genes":len(selected),
        "selected_all_isoforms":len(all_selected_sequences),"selected_all_isoforms_sha256":sha256(selected_fasta),
        "original_training_isoforms":len(train_accessions),"original_training_fasta_sha256":sha256(train_fasta),
        "candidate_order_sha256":sha256(out / "candidate_order.jsonl"),"selection_ledger_sha256":sha256(out / "selection_ledger.jsonl"),
        "homology_edges_sha256":sha256(out / "homology_edges.jsonl"),"source_check":checks,
        "mmseqs_version":subprocess.check_output([str(tool),"version"],text=True).strip(),"mmseqs_sha256":sha256(tool),
        "searches":commands,"prepare_script_sha256":sha256(Path(__file__)),"rejected_candidates":rejected}
    write_json(out / "manifest.json",manifest)
    print(json.dumps({"status":manifest["status"],"targets":len(selected),"checks":checks}),flush=True)


if __name__=="__main__":main()
