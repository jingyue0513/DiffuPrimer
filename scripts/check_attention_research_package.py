#!/usr/bin/env python3
"""Check archive hashes and run the public example from actual extracted files."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive",default="paper/DiffuPrimer_attention_research_materials.zip")
    p.add_argument("--python",default=str(ROOT/".venv_attention/bin/python"))
    p.add_argument("--qprimer-python",default=str(ROOT/".venv/bin/python"))
    args=p.parse_args()
    archive=ROOT/args.archive
    checked=0
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError("Archive CRC failure")
        manifest=json.loads(z.read("DiffuPrimer/research_materials/package_manifest.json"))
        for row in manifest["files"]:
            data=z.read("DiffuPrimer/"+row["path"])
            if len(data)!=row["bytes"] or hashlib.sha256(data).hexdigest()!=row["sha256"]:
                raise ValueError("Archive member hash failure: "+row["path"])
            checked+=1
        with tempfile.TemporaryDirectory(prefix="diffuprimer-zip-check-") as tmp:
            prefix=("src/","scripts/","configs/","experiments/hybrid_attention/src/",
                    "experiments/hybrid_attention/results/checkpoints/hybrid_attention_diffusion_seed20261013.pt",
                    "experiments/attention/baselines/run_qprimer.py",
                    "experiments/attention/baselines/vendor/qPrimer/","examples/")
            for name in z.namelist():
                rel=name.removeprefix("DiffuPrimer/")
                if not rel.startswith(prefix):continue
                path=Path(tmp)/rel
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_bytes(z.read(name))
            env=dict(os.environ);env.pop("PYTHONPATH",None);env["CUDA_VISIBLE_DEVICES"]=""
            command=[args.python,"scripts/design_attention.py","--fasta","examples/example.fasta","--output","result.json","--budget","2","--device","cpu"]
            subprocess.run(command,cwd=tmp,env=env,check=True,capture_output=True,text=True)
            result=json.loads((Path(tmp)/"result.json").read_text())
            valid=sum(bool(c.get("structural_valid") and c.get("assessment",{}).get("computational_feasible")) for c in result["candidates"])
            if len(result["candidates"])!=2 or valid<1:raise ValueError("Extracted design CLI returned no feasible pair")
            # load_upstream verifies both commit provenance and genuine design.py SHA,
            # including the portable source snapshot when .git is intentionally absent.
            probe="import runpy; m=runpy.run_path('experiments/attention/baselines/run_qprimer.py'); m['load_upstream']()"
            subprocess.run([args.qprimer_python,"-c",probe],cwd=tmp,env=env,check=True,capture_output=True,text=True)
    record={"actual_archive_extracted":True,"archive_members_sha256_verified":checked,
            "cpu_design_valid_pairs":valid,"pinned_qprimer_without_git_loads":True,
            "archive_sha256":hashlib.sha256(archive.read_bytes()).hexdigest()}
    (ROOT/"paper/attention_research_package_extract_check.json").write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps(record))


if __name__=="__main__":main()
