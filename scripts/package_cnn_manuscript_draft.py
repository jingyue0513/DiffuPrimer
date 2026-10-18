#!/usr/bin/env python3
"""Package the checked, editable research draft without claiming submission."""
import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'paper/bioinformatics_cnn_draft'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    check=json.loads((OUT/'verification.json').read_text())
    assert check['passed'] and check['main_pages']<=7
    assert sha(OUT/'manuscript.pdf')==check['manuscript_pdf_sha256']
    assert sha(OUT/'supplementary.pdf')==check['supplementary_pdf_sha256']
    suffixes={'.tex','.pdf','.json','.bib','.cls','.bst','.md','.txt','.png','.tiff'}
    files=[p for p in OUT.rglob('*') if p.is_file() and p.suffix in suffixes
           and ('previews' not in p.parts or p.name.endswith('_contact_sheet.png'))
           and p.name not in {'main_05_mechanism.pdf','preview_first_page.png'}]
    files.extend([ROOT/'scripts/build_cnn_manuscript_materials.py',ROOT/'scripts/check_cnn_manuscript_materials.py',
                  ROOT/'scripts/package_cnn_manuscript_draft.py',ROOT/'docs/bioinformatics_submission_requirements.md',
                  ROOT/'docs/predac_manuscript_style_audit.json',ROOT/'docs/predac_manuscript_style_audit.md',
                  ROOT/'paper/figures/architecture_method/DiffuPrimer_method.pptx',
                  ROOT/'paper/figures/architecture_method/DEJAVU_FONT_LICENSE.txt'])
    records=[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]
    bundle=ROOT/'paper/DiffuPrimer_English_research_draft.zip'
    manifest={'at':datetime.now(timezone.utc).isoformat(),'document_status':'AI-assisted research draft for author review',
              'journal_policy_compliance_asserted':False,'main_pages':check['main_pages'],'supplementary_pages':check['supplementary_pages'],
              'references':52,'files':records}
    with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(files):z.write(p,'DiffuPrimer_draft/'+str(p.relative_to(ROOT)))
        z.writestr('DiffuPrimer_draft/bundle_manifest.json',json.dumps(manifest,indent=2)+'\n')
    with zipfile.ZipFile(bundle) as z:
        assert z.testzip() is None
        for record in records:
            payload=z.read('DiffuPrimer_draft/'+record['path'])
            assert len(payload)==record['bytes'] and hashlib.sha256(payload).hexdigest()==record['sha256']
    receipt={'at':manifest['at'],'archive':str(bundle.relative_to(ROOT)),'bytes':bundle.stat().st_size,'sha256':sha(bundle),
             'files_checked':len(files),'all_entry_hashes_passed':True,'crc_passed':True,
             'main_pages':check['main_pages'],'supplementary_pages':check['supplementary_pages'],'references':52}
    (ROOT/'paper/cnn_manuscript_draft_bundle_check.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))


if __name__=='__main__':main()
