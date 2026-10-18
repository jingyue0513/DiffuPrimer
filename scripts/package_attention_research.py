#!/usr/bin/env python3
"""Build a local, portable archive of completed computational research materials.

Full candidates are byte-preserving copies. Reduced screen exports are explicitly
named exports and retain the SHA-256 of the original complete witness records.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = Path("experiments/attention/results/study_summary.json")
EXTERNAL = Path("experiments/attention_external/results/study_summary.json")
CHECKPOINTS = (
    "results/checkpoints/diffusion_weak_partial_pilot_seed20261013.pt",
    "results/checkpoints/direct_weak_partial_pilot_seed20261013.pt",
    "experiments/attention/results/checkpoints/attention_diffusion_seed20261013.pt",
    "experiments/attention/results/checkpoints/attention_direct_seed20261013.pt",
    "experiments/hybrid_attention/results/checkpoints/hybrid_attention_diffusion_seed20261013.pt",
    "experiments/hybrid_attention/results/checkpoints/hybrid_attention_direct_seed20261013.pt",
)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def relative(path):
    p = Path(path)
    return p.resolve().relative_to(ROOT) if p.is_absolute() else p


def dump(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+"\n")


def source_files():
    selected = set()
    for folder in ("src","scripts","tests","configs"):
        selected.update(p for p in (ROOT/folder).rglob("*") if p.is_file() and p.suffix in {".py",".cjs",".json",".yaml",".yml"} and "__pycache__" not in p.parts)
    for p in (ROOT/"experiments").rglob("*"):
        if not p.is_file() or "__pycache__" in p.parts or "vendor" in p.parts or "script_history" in p.parts:
            continue
        if "results" not in p.parts and p.suffix in {".py",".json",".txt",".md",".toml"}:
            selected.add(p)
    for name in ("README.md","pyproject.toml","requirements.lock",".gitignore"):
        selected.add(ROOT/name)
    for p in (ROOT/"docs").glob("*.md"):
        selected.add(p)
    for p in (ROOT/"data/splits").glob("*.json*"):
        selected.add(p)
    for name in ("dataset_partial_pilot.jsonl","dataset_partial_pilot.manifest.json","targets_pilot.jsonl",
                 "attention_external_targets.jsonl","pilot_gene_selection.jsonl","pilot_gene_selection.manifest.json"):
        selected.add(ROOT/"data/processed"/name)
    for name in ("results/reference_verification.json","results/spliced_reference_verification.json",
                 "results/checks/span_cache/check.json","results/checks/bowtie_exhaustive/check.json",
                 "results/checks/sqlite_specificity/check.json","tools/manifest.json","tools/bowtie_manifest.json",
                 "paper/author_metadata.json","experiments/final_followup/redesign_tasks.jsonl",
                 "experiments/attention_external/cohort.csv",
                 "experiments/attention_external/candidate_order.jsonl","experiments/attention_external/selection_ledger.jsonl",
                 "experiments/attention_external/gene_selection.jsonl","experiments/attention_external/external_split.jsonl",
                 "experiments/attention_external/homology_edges.jsonl","experiments/attention_external/selected_all_isoforms.fasta",
                 "data/processed/pilot_all_isoforms.fasta"):
        selected.add(ROOT/name)
    selected.update(ROOT/x for x in CHECKPOINTS)
    # Exact upstream source and existing license texts, without git history,
    # generated outputs, native indices, or an invented license for original code.
    vendor = ROOT/"experiments/attention/baselines/vendor"
    for child in ("qPrimer","PrimerCast","sam2pairwise"):
        for p in (vendor/child).rglob("*"):
            if not p.is_file() or ".git" in p.parts or "__pycache__" in p.parts or "tests" in p.parts:
                continue
            if p.suffix in {".py",".toml",".yml",".txt",".pth",".joblib",".cc",".h",".json",".md"} or p.name in {"LICENSE","Makefile","README"}:
                selected.add(p)
    # Native Linux x86_64 tools make the recorded full-reference pipeline replayable.
    for base,names in {
        "tools/bowtie-1.3.1-linux-x86_64":("bowtie","bowtie-build","bowtie-align-s","bowtie-align-l","bowtie-build-s","bowtie-build-l","LICENSE"),
        "tools/ncbi-blast-2.17.0+":("bin/blastn","bin/makeblastdb","LICENSE","README"),
        "tools/mmseqs":("bin/mmseqs","LICENSE.md","README.md"),
    }.items():
        selected.update(ROOT/base/name for name in names)
    # Original source documentation and official authoring-template files.
    for directory in ("paper/bioinformatics_template",):
        for p in (ROOT/directory).rglob("*"):
            if p.is_file() and p.suffix in {".cls",".bst",".bib",".tex",".json",".md",".pdf",".xlsx"}:
                selected.add(p)
    return {p for p in selected if p.is_file()}


def validate_final():
    documents={"main":read_json(ROOT/MAIN),"external":read_json(ROOT/EXTERNAL)}
    for label,doc in documents.items():
        if doc.get("status")!="complete" or not doc.get("additional_experiments_complete") or doc.get("allow_partial") or not doc.get("inputs"):
            raise ValueError(f"{label} analysis is not a completed final source set")
        if doc.get("skipped_sources") or doc.get("pending_planned_runs"):
            raise ValueError(f"{label} analysis still has pending or skipped sources")
    final_protocol=ROOT/"experiments/attention/final_protocol.json"
    if documents["main"].get("protocol_sha256")!=digest(final_protocol):
        raise ValueError("Main analysis has not been rebuilt with the final architecture protocol")
    for label,summary in (("main",MAIN),("external",EXTERNAL)):
        provenance=ROOT/summary.parent/"figures/provenance.json"
        record=read_json(provenance)
        if record.get("summary_sha256")!=digest(ROOT/summary) or record.get("status")!="complete":
            raise ValueError(f"{label} publication plots do not match the completed final summary")
        for figure in record.get("generated",[]):
            for asset in figure.get("files",[]):
                if digest(ROOT/relative(asset["path"]))!=asset["sha256"]:
                    raise ValueError(f"Publication figure bytes changed: {asset['path']}")
    page=ROOT/"reports/attention_research"
    page_manifest=read_json(page/"build_manifest.json")
    if page_manifest.get("presentation_status")!="completed_results_review" or page_manifest.get("summary_sha256")!=digest(ROOT/MAIN) or page_manifest.get("external_summary_sha256")!=digest(ROOT/EXTERNAL):
        raise ValueError("Final research webpage has not been rebuilt from both completed analyses")
    if page_manifest.get("page_sha256")!=digest(page/"index.html"):
        raise ValueError("Research webpage bytes do not match its build manifest")
    browser=read_json(page/"browser_check.json")
    if not browser.get("passed") or browser.get("page_sha256")!=digest(page/"index.html") or browser.get("summary_sha256")!=digest(ROOT/MAIN):
        raise ValueError("Final research webpage has not passed browser checks for its current bytes")
    worksheet=ROOT/"paper/bioinformatics_attention_materials"
    worksheet_manifest=read_json(worksheet/"source_manifest.json")
    if worksheet_manifest.get("mode")!="completed-source author result material" or worksheet_manifest.get("main_summary_sha256")!=digest(ROOT/MAIN) or worksheet_manifest.get("external_summary_sha256")!=digest(ROOT/EXTERNAL):
        raise ValueError("Author result worksheet has not been rebuilt from both completed analyses")
    layout=read_json(worksheet/"layout_report.json")
    if not layout.get("within_target") or not layout.get("numerical_result_instances") or layout.get("pdf_sha256")!=digest(worksheet/"author_results.pdf"):
        raise ValueError("Author result worksheet PDF does not match the final measured layout")
    for figure in worksheet_manifest.get("figures",[]):
        if digest(ROOT/relative(figure["source"]))!=figure["sha256"]:
            raise ValueError(f"Author worksheet figure source changed: {figure['source']}")
        if digest(worksheet/"figures"/(figure["name"]+".pdf"))!=figure["sha256"]:
            raise ValueError(f"Author worksheet copied figure changed: {figure['name']}")
    delivery=read_json(ROOT/"results/attention_delivery_check.json")
    if not delivery.get("computational_research_complete") or not delivery.get("browser_checked"):
        raise ValueError("Final research delivery has not passed artifact verification")
    for artifact in delivery["artifacts"]:
        if digest(ROOT/artifact["path"])!=artifact["sha256"]:
            raise ValueError(f"Final delivered artifact changed: {artifact['path']}")
    return documents


def smoke_check():
    """Exercise the actual distribution CLI from an isolated copied directory."""
    destination=ROOT/"results/package_smoke_check"
    destination.mkdir(parents=True,exist_ok=True)
    targets=[json.loads(line) for line in (ROOT/"data/processed/targets_pilot.jsonl").read_text().splitlines()]
    target=next(r for r in targets if r["gene_id"]=="16716")
    fasta=f">{target['transcript_accession_version']} example_first_512_bases\n{target['template_sequence'][:512]}\n"
    (destination/"example.fasta").write_text(fasta)
    with tempfile.TemporaryDirectory(prefix="diffuprimer-portable-") as tmp:
        copied=Path(tmp)
        for directory in ("src","scripts","configs","experiments/hybrid_attention/src"):
            shutil.copytree(ROOT/directory,copied/directory,ignore=shutil.ignore_patterns("__pycache__"))
        checkpoint=Path(CHECKPOINTS[-2])
        (copied/checkpoint.parent).mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/checkpoint,copied/checkpoint)
        shutil.copy2(destination/"example.fasta",copied/"example.fasta")
        python=ROOT/".venv_attention/bin/python"
        if not python.exists():python=Path(sys.executable)
        command=[str(python),"scripts/design_attention.py","--fasta","example.fasta","--output","example_result.json","--budget","2","--device","cpu"]
        env=dict(os.environ);env.pop("PYTHONPATH",None);env["CUDA_VISIBLE_DEVICES"]=""
        with (destination/"stdout.txt").open("w") as log:
            subprocess.run(command,cwd=copied,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        result=read_json(copied/"example_result.json")
        if len(result["candidates"])!=2 or not any(c.get("structural_valid") and c.get("assessment",{}).get("computational_feasible") for c in result["candidates"]):
            raise ValueError("Portable CLI example returned no valid primer-pair proposal")
        shutil.copy2(copied/"example_result.json",destination/"example_result.json")
    record={"passed":True,"command":command,"execution_directory":"isolated temporary package-style copy",
        "python":str(python),"design_script_sha256":digest(ROOT/"scripts/design_attention.py"),
        "checkpoint_sha256":digest(ROOT/CHECKPOINTS[-2]),"example_sha256":digest(destination/"example.fasta"),
        "result_sha256":digest(destination/"example_result.json"),"budget":2,"device":"cpu"}
    dump(destination/"check.json",record)
    return record


def package(documents,smoke,archive):
    files=source_files()
    for directory in ("reports/attention_research","paper/bioinformatics_attention_materials"):
        files.update(p for p in (ROOT/directory).rglob("*") if p.is_file() and p.suffix in {".html",".pdf",".svg",".png",".pptx",".json",".csv",".tex",".bib",".cls",".bst",".md",".txt",".otf"} and not p.name.startswith("preview_"))
    for name in ("experiments/attention/results/status.json","experiments/attention/results/final_verification.log",
                 "experiments/attention_external/results/final_verification.log",
                 "experiments/attention/results/model_check.json","experiments/attention/baselines/qprimer_generation_summary.json",
                 "experiments/attention/baselines/primercast_generation_summary.json",
                 "experiments/attention_external/results/baseline_generation_summary.json",
                 "results/current_research_status.json","results/attention_delivery_check.json"):
        path=ROOT/name
        if path.is_file():files.add(path)
    for summary in (MAIN,EXTERNAL):
        files.update(p for p in (ROOT/summary.parent).glob("study_summary.*") if p.is_file())
        files.update(p for p in (ROOT/summary.parent/"figures").glob("*") if p.is_file() and p.suffix in {".pdf",".svg",".png",".json"})
    for directory in ("paper/figures/hybrid_attention_method","paper/figures/exon_interpretation"):
        files.update(p for p in (ROOT/directory).glob("*") if p.is_file() and p.suffix in {".pptx",".pdf",".svg",".png",".tiff",".json",".md"})
    for directory in ("experiments/attention/results/exon_interpretation","experiments/attention/results/pickr_literature"):
        files.update(p for p in (ROOT/directory).glob("*") if p.is_file() and p.suffix in {".json",".jsonl",".csv"})
    for directory in ("experiments/attention/results","experiments/hybrid_attention/results"):
        files.update(p for p in (ROOT/directory).glob("*training.json") if p.is_file())
        selected=ROOT/directory/"validation_quality_selection.json"
        if selected.exists():files.add(selected)
    with tempfile.TemporaryDirectory(prefix="diffuprimer-materials-") as tmp:
        generated=Path(tmp)
        (generated/"research_materials").mkdir()
        dump(generated/"experiments/attention/baselines/vendor/qPrimer/source_snapshot.json",{
            "repository":"https://github.com/swu1019lab/qPrimer",
            "commit":"34f047e5d77ec0353846a9a18b0b94f74013db3b",
            "design_sha256":digest(ROOT/"experiments/attention/baselines/vendor/qPrimer/qPrimer/design.py"),
            "acquisition":"Original pinned git source copied without git history for the portable archive"})
        exports=[];replay_sources=[];seen=set()
        for label,doc in documents.items():
            for source in [*doc["inputs"],*doc.get("validation_inputs",[])]:
                directory=ROOT/relative(source["directory"])
                rel=directory.relative_to(ROOT)
                if str(rel) in seen:continue
                seen.add(str(rel))
                required=["candidates.jsonl","per_gene.jsonl","manifest.json","completion.json","specificity_bowtie_status.json"]
                for name in required:
                    path=directory/name
                    if not path.exists():raise ValueError(f"Missing completed source record {path}")
                    expected=source.get(name+"_sha256")
                    if expected and digest(path)!=expected:raise ValueError(f"Pinned source bytes changed: {path}")
                    files.add(path)
                screen_path=directory/"specificity_bowtie.jsonl"
                actual_screen_sha=digest(screen_path)
                if source.get("specificity_bowtie.jsonl_sha256") and actual_screen_sha!=source["specificity_bowtie.jsonl_sha256"]:
                    raise ValueError(f"Pinned full screen bytes changed: {screen_path}")
                screens={}
                with screen_path.open() as f:
                    for line in f:
                        row=json.loads(line)
                        keep={key:row[key] for key in ("candidate_id","status","gene_specificity_screen_pass","genomic_DNA_product_detected","genome_status","scope") if key in row}
                        keep["product_counts"]={key:len(row.get(key,[])) for key in ("on_target_products","off_target_products","unresolved_products","genomic_potential_amplicons")}
                        screens[row["candidate_id"]]=keep
                export_root=generated/"research_materials/compact_runs"/rel
                export_root.mkdir(parents=True,exist_ok=True)
                candidate_path=directory/"candidates.jsonl"
                joined_path=export_root/"candidate_screen.export.jsonl"
                outcomes_path=export_root/"screen_outcomes.export.jsonl"
                with candidate_path.open() as inputs,joined_path.open("w") as joined:
                    count=0
                    for line in inputs:
                        candidate=json.loads(line);cid=candidate["candidate_id"]
                        if cid not in screens:raise ValueError(f"Unscreened candidate in completed source: {cid}")
                        joined.write(json.dumps({"candidate":candidate,"screen_outcome":screens[cid]},separators=(",",":"))+"\n");count+=1
                with outcomes_path.open("w") as f:
                    for row in screens.values():f.write(json.dumps(row,separators=(",",":"))+"\n")
                export_record={"source_directory":str(directory),"portable_source_directory":str(rel),"analysis_cohort":label,
                    "method":source.get("method"),"source_family":source.get("source_family"),"rows":count,
                    "original_candidates_sha256":digest(candidate_path),"original_complete_screen_sha256":actual_screen_sha,
                    "export_schema":"candidate_screen.export/v1",
                    "export_description":"Candidate fields plus screen flags and product counts; complete product-coordinate witnesses are omitted. This is a derived export, not the original full screening record.",
                    "candidate_export_sha256":digest(joined_path),"outcome_export_sha256":digest(outcomes_path),
                    "original_candidate_record_copied_byte_for_byte":True,"complete_screen_witness_record_included":False}
                dump(export_root/"export_provenance.json",export_record)
                exports.append(export_record)
                replay_sources.append({"directory":str(rel),"cohort":label,"method":source.get("method"),"source_family":source.get("source_family"),"validation":source not in doc["inputs"]})
        # Index manifests are evidence, not a claim that absent indices are installed.
        for base in (ROOT/"data/processed/bowtiedb",ROOT/"data/processed/blastdb"):
            for p in base.glob("*.manifest.json"):
                copied=generated/"research_materials/reference_manifests"/base.name/p.name
                copied.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,copied)
        dump(generated/"research_materials/replay_sources.json",replay_sources)
        dump(generated/"research_materials/export_index.json",exports)
        example=generated/"examples"
        example.mkdir()
        for name in ("example.fasta","example_result.json","check.json"):
            shutil.copy2(ROOT/"results/package_smoke_check"/name,example/name)
        shutil.copy2(ROOT/"docs/attention_research_materials_README.md",generated/"RESEARCH_MATERIALS_README.md")
        manifest={"package_kind":"completed computational research materials; local archive",
            "formal_paper_completion_asserted":False,"public_upload_performed":False,
            "main_summary_sha256":digest(ROOT/MAIN),"external_summary_sha256":digest(ROOT/EXTERNAL),
            "main_status":documents["main"]["status"],"external_status":documents["external"]["status"],
            "original_training_dataset_sha256":digest(ROOT/"data/processed/dataset_partial_pilot.jsonl"),
            "checkpoint_files":[{"path":p,"sha256":digest(ROOT/p)} for p in CHECKPOINTS],
            "smoke_check":smoke,"exported_source_runs":len(exports),
            "excluded":["raw full-reference SAM alignments","alignment SQLite cache","native/vendor indices and temporary files","large raw RNA/genome references","full product-coordinate screening witnesses"],
            "licensing":"Original project license is not assigned by this builder. Existing upstream license texts and copyright headers are preserved.",
            "files":[]}
        members={str(p.relative_to(ROOT)):p for p in files}
        members.update({str(p.relative_to(generated)):p for p in generated.rglob("*") if p.is_file()})
        for name,path in sorted(members.items()):
            manifest["files"].append({"path":name,"bytes":path.stat().st_size,"sha256":digest(path)})
        dump(generated/"research_materials/package_manifest.json",manifest)
        members["research_materials/package_manifest.json"]=generated/"research_materials/package_manifest.json"
        archive.parent.mkdir(parents=True,exist_ok=True)
        temporary=archive.with_suffix(archive.suffix+".tmp")
        with zipfile.ZipFile(temporary,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as zipout:
            for name,path in sorted(members.items()):zipout.write(path,"DiffuPrimer/"+name)
        temporary.replace(archive)
        with zipfile.ZipFile(archive) as zin:
            if zin.testzip() is not None:raise ValueError("Archive CRC verification failed")
        receipt={"archive":str(archive.relative_to(ROOT)),"bytes":archive.stat().st_size,"sha256":digest(archive),
            "entries":len(members),"source_runs":len(exports),"main_status":"complete","external_status":"complete",
            "archive_crc_verified":True,"cpu_example_passed":True,"formal_paper_completion_asserted":False}
        dump(ROOT/"paper/attention_research_package_manifest.json",receipt)
        print(json.dumps(receipt,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",default="paper/DiffuPrimer_attention_research_materials.zip")
    p.add_argument("--smoke-test-only",action="store_true",help="Test portable custom-sequence design without creating a final archive")
    args=p.parse_args()
    if args.smoke_test_only:
        print(json.dumps(smoke_check(),indent=2));return
    documents=validate_final()
    smoke=smoke_check()
    package(documents,smoke,ROOT/args.output)


if __name__=="__main__":main()
