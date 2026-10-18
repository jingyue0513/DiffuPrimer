#!/usr/bin/env python3
"""Acquire fixed upstream source and pretrained weights without git history."""
import concurrent.futures
import hashlib
import json
from pathlib import Path

import requests

COMMIT = "50d357b43a727196b5ae2e1c2e7d041519243673"
REPO = "broadinstitute/PrimerCast"
DEST = Path(__file__).resolve().parent / "vendor/PrimerCast"


def main():
    r = requests.get(f"https://api.github.com/repos/{REPO}/git/trees/{COMMIT}?recursive=1",timeout=60)
    r.raise_for_status()
    metadata = r.json()
    wanted = [x["path"] for x in metadata["tree"] if x["type"]=="blob" and
        (x["path"].startswith("src/primercast/") or x["path"] in
         ["README.md","LICENSE","pyproject.toml","environment.yml","workflows/params.txt.template","docs/primercast_cli.md"])]
    def fetch(path):
        response = requests.get(f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/{path}",timeout=60)
        response.raise_for_status()
        dest = DEST / path
        dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_bytes(response.content)
        return {"path":path,"bytes":len(response.content),"sha256":hashlib.sha256(response.content).hexdigest()}
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        records = list(executor.map(fetch,wanted))
    (DEST / "source_snapshot.json").write_text(json.dumps({"repository":"https://github.com/"+REPO,
        "commit":COMMIT,"files":records},indent=2))
    print(json.dumps({"commit":COMMIT,"files":len(records)}))


if __name__=="__main__":
    main()
