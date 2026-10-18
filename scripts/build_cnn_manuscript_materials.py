#!/usr/bin/env python3
"""Derive manuscript numbers and supplemental tables from frozen study results."""
import hashlib
import io
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'paper/bioinformatics_cnn_draft'
COUNT = 'unique_RNA_pass_pairs_without_genomic_flag'
RNA = 'unique_RNA_pass_pairs'
FRACTION = 'RNA_no_genome_fraction'
DIVERSITY = 'RNA_pass_amplicon_distance'
THERMAL = 'best_RNA_pass_temperature_penalty'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source_paths = [ROOT / 'experiments/attention/results/study_summary.json',
                    ROOT / 'experiments/attention_external/results/study_summary.json']
    main_data, extra = [json.loads(p.read_text()) for p in source_paths]
    assert all(s['status'] == 'complete' and s['total_test_genes'] == 75 for s in (main_data, extra))
    records, macros = [], []

    def group(data, method, context='amp80-200', variant='full', budget=20):
        return next(g for g in data['summaries'] if (g['method'], g['context'], g['variant'], g['budget']) == (method, context, variant, budget))

    def value(name, data, method, endpoint=COUNT, context='amp80-200', variant='full', budget=20, scale=1, digits=2):
        g = group(data, method, context, variant, budget)
        v = g[endpoint]
        assert v is not None
        macros.append('\\newcommand{\\' + name + '}{' + f'{v*scale:.{digits}f}' + '}')
        records.append(dict(macro=name, cohort='original' if data is main_data else 'additional', context=context,
                            method=method, variant=variant, budget=budget, endpoint=endpoint, value=v,
                            scale=scale, digits=digits, n=g[endpoint + '_n'], ci95=g[endpoint + '_ci95'], run=g['run']))
        return v

    def contrast(name, data, right, endpoint=COUNT, context='amp80-200', scale=1, digits=2):
        c = next(c for c in data['paired_comparisons'] if
                 (c['context'],c['left_method'],c['right_method'],c['left_variant'],c['right_variant'],c['budget'],c['endpoint']) ==
                 (context,'cnn_diffusion',right,'full','full',20,endpoint))
        for suffix, v in [('Delta',c['mean_difference']),('Lo',c['ci95'][0]),('Hi',c['ci95'][1])]:
            macros.append('\\newcommand{\\' + name + suffix + '}{' + f'{v*scale:.{digits}f}' + '}')
        macros.append('\\newcommand{\\' + name + 'N}{' + str(c['n']) + '}')
        records.append(dict(macro_prefix=name, cohort='original' if data is main_data else 'additional',
                            contrast=c, scale=scale, digits=digits))

    for label, data in [('Main',main_data),('Extra',extra)]:
        for method, label_m in [('cnn_diffusion','CNN'),('cnn_direct','Direct'),('primer3','Primer'),('qprimer','QPrimer'),
                                ('attention_rope_diffusion','Attention'),('attention_hybrid_diffusion','Hybrid')]:
            value(label+label_m+'Ret',data,method)
        for method,label_m in [('cnn_diffusion','CNN'),('cnn_direct','Direct'),('primer3','Primer'),('qprimer','QPrimer')]:
            value(label+label_m+'RNA',data,method,RNA)
            value(label+label_m+'Fraction',data,method,FRACTION,scale=100,digits=1)
        for method,label_m in [('cnn_diffusion','CNN'),('cnn_direct','Direct')]:
            value(label+label_m+'Diversity',data,method,DIVERSITY,digits=3)
            value(label+label_m+'Thermal',data,method,THERMAL,digits=3)
        contrast(label+'Primer',data,'primer3')
        contrast(label+'PrimerFraction',data,'primer3',FRACTION,scale=100)
        contrast(label+'Direct',data,'cnn_direct')
        contrast(label+'DirectDiversity',data,'cnn_direct',DIVERSITY,digits=3)
        contrast(label+'DirectThermal',data,'cnn_direct',THERMAL,digits=3)
        a,b=group(data,'cnn_diffusion')[COUNT],group(data,'primer3')[COUNT]
        gain=(a/b-1)*100
        macros.append('\\newcommand{\\'+label+'PrimerGain}{'+f'{gain:.1f}'+'}')
        records.append(dict(macro=label+'PrimerGain',derivation='100*(CNN retained mean / Primer3 retained mean - 1)',value=gain))
    for context,tag in [('amp80-120','Short'),('amp121-160','Middle'),('amp161-200','Long'),
                        ('amp80-200_fix_forward','CompleteF'),('amp80-200_fix_reverse','CompleteR'),
                        ('redesign_forward','RedesignF'),('redesign_reverse','RedesignR')]:
        for m,t in [('cnn_diffusion','CNN'),('cnn_direct','Direct'),('primer3','Primer'),('qprimer','QPrimer')]:
            value(tag+t+'Ret',main_data,m,context=context)
        value(tag+'CNNRNA',main_data,'cnn_diffusion',RNA,context=context)
    for variant,tag in [('geometry_only','Geometry'),('single_oligo','Oligo'),('steps1','OneStep'),('steps4','FourStep'),('steps8','EightStep')]:
        value(tag+'Ret',main_data,'cnn_diffusion',variant=variant)
        value(tag+'RNA',main_data,'cnn_diffusion',RNA,variant=variant)
    (OUT/'result_values.tex').write_text('% Generated only from completed frozen summaries.\n'+'\n'.join(macros)+'\n')

    template=ROOT/'paper/bioinformatics_template/oup-authoring-template'
    for fn in ['oup-authoring-template.cls','oup-abbrvnat.bst']:
        shutil.copyfile(template/fn,OUT/fn)
    figures={
      '01_method':'paper/figures/architecture_method/DiffuPrimer_method.pdf',
      '02_quality':'paper/figures/cnn_publication_results/fig2_candidate_quality.pdf',
      '03_flexibility':'paper/figures/cnn_publication_results/fig3_flexible_design.pdf',
      '04_biology':'paper/figures/cnn_biological_analysis/fig4_biological_specificity.pdf',
      '05_mechanism':'paper/figures/cnn_publication_results/fig5_mechanism.pdf',
    }
    (OUT/'figures').mkdir(exist_ok=True)
    for name,source in figures.items():shutil.copyfile(ROOT/source,OUT/'figures'/f'{name}.pdf')
    # Main-article subsets retain vector data and native text size. Full panels
    # remain in the supplementary figures and the complete analysis webpage.
    crops={'02_quality':(.53,1), '03_flexibility':(.36,1), '04_biology':(0,.655)}
    crop_records=[]
    for name,(lo,hi) in crops.items():
        original=OUT/'figures'/f'{name}.pdf';pg=PdfReader(original).pages[0]
        w,h=float(pg.mediabox.width),float(pg.mediabox.height);bottom,top=lo*h,hi*h
        pg.cropbox=RectangleObject([0,bottom,w,top]);writer=PdfWriter()
        target=writer.add_blank_page(width=w,height=top-bottom)
        target.merge_transformed_page(pg,Transformation().translate(0,-bottom))
        if name=='04_biology':
            fig=plt.figure(figsize=(w/72,(top-bottom)/72));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,w),ylim=(0,top-bottom));ax.axis('off')
            for old,new,y in [('B','A',310.5354330709),('C','B',222.7322834646),('D','C',114.094488189)]:
                ax.add_patch(Rectangle((16.5,y-2),10,13,facecolor='white',edgecolor='none'))
                ax.text(17.6598425197,y,new,fontsize=10,fontweight='bold',color='#193D50',va='baseline',ha='left')
            buffer=io.BytesIO();fig.savefig(buffer,format='pdf',transparent=True);plt.close(fig);buffer.seek(0)
            target.merge_page(PdfReader(buffer).pages[0])
        output=OUT/'figures'/f'main_{name}.pdf'
        with output.open('wb') as f:writer.write(f)
        crop_records.append({'source':f'figures/{name}.pdf','source_sha256':sha(original),'main_output':f'figures/main_{name}.pdf','main_sha256':sha(output),'crop_y_fraction':[lo,hi],'width_unchanged':True,'biological_panel_relabelling':'B/C/D to A/B/C' if name=='04_biology' else None})
    labels={'cnn_diffusion':'DiffuPrimer (CNN diffusion)','cnn_direct':'CNN direct','primer3':'Primer3','qprimer':'qPrimer design',
            'primercast':'PrimerCast','uniform_pair_sampler':'Uniform pair','rule_sampler':'Forward-first',
            'attention_rope_diffusion':'RoPE diffusion','attention_rope_direct':'RoPE direct',
            'attention_hybrid_diffusion':'CNN + attention diffusion','attention_hybrid_direct':'CNN + attention direct'}
    table_methods=[m for m in labels if any(g['method']==m and g['context']=='amp80-200' and g['variant']=='full' for g in main_data['summaries'])]
    def tex_table(name, caption, headings, rows, columns):
        return '\\begin{table*}[t]\\centering\\small\\setlength{\\tabcolsep}{3pt}\n\\caption{'+caption+'}\\label{tab:'+name+'}\n\\begin{tabular}{@{}'+columns+'@{}}\\toprule\n'+' & '.join(headings)+' \\\\\\midrule\n'+'\n'.join(' & '.join(row)+' \\\\' for row in rows)+'\n\\bottomrule\\end{tabular}\\end{table*}\n'
    rows=[]
    for m in table_methods:
        row=[labels[m]]
        for data in [main_data,extra]:
            try:g=group(data,m)
            except StopIteration:row.extend(['--']*3);continue
            row.extend([f'{g[RNA]:.2f}',f'{g[COUNT]:.2f}',f'{100*g[FRACTION]:.1f} ({g[FRACTION+"_n"]})' if g[FRACTION] is not None else '--'])
        rows.append(row)
    tables=tex_table('allmethods','Complete-design results at $K=20$, with the original and additional cohorts reported separately. RNA and retained are mean unique pairs per gene ($N=75$ in each cohort). Fraction is the mean retained/RNA percentage with defined-gene count in parentheses. A dash denotes an experiment not performed.',
                    ['Method','Original RNA','Retained','Fraction (\\%)','Additional RNA','Retained','Fraction (\\%)'],rows,'lrrrrrr')
    rows=[]
    for context,tag in [('amp80-120','80--120 bp'),('amp121-160','121--160 bp'),('amp161-200','161--200 bp'),('amp80-200_fix_forward','Completion: retain F'),('amp80-200_fix_reverse','Completion: retain R'),('redesign_forward','Redesign: retain F'),('redesign_reverse','Redesign: retain R')]:
        rows.append([tag]+[f'{group(main_data,m,context)[COUNT]:.2f}' for m in ['cnn_diffusion','cnn_direct','primer3','qprimer']])
    tables+=tex_table('tasks','Retained-pair counts for length control, completion and changed-length redesign. All 75 original test genes are included; $K=20$.',
                      ['Task','DiffuPrimer','CNN direct','Primer3','qPrimer design'],rows,'lrrrr')
    rows=[]
    for variant,tag in [('raw','Raw independent labels'),('geometry_only','Geometry support'),('single_oligo','Single-oligo guidance'),('full','Full pair guidance'),('steps1','1 step'),('steps4','4 steps'),('steps8','8 steps')]:
        g=group(main_data,'cnn_diffusion',variant=variant)
        rows.append([tag]+[f'{g[k]:.3f}' if g[k] is not None else '--' for k in [RNA,COUNT,DIVERSITY,THERMAL]])
    tables+=tex_table('components','Same-checkpoint guidance and reverse-step experiments. Diversity is the mean RNA-specific amplicon Jaccard distance and thermal is the best RNA-specific temperature penalty. Endpoint-specific denominators are recorded in the supplied summary.',
                      ['Configuration','RNA','Retained','Diversity','Thermal'],rows,'lrrrr')
    rows=[]
    for data,tag in [(main_data,'Original'),(extra,'Additional')]:
        for endpoint,title,scale in [(COUNT,'Retained count',1),(FRACTION,'Retained fraction (pp)',100),(DIVERSITY,'Jaccard distance',1),(THERMAL,'Temperature penalty',1)]:
            for right in ['cnn_direct','primer3','qprimer']:
                c=next(c for c in data['paired_comparisons'] if (c['context'],c['left_method'],c['right_method'],c['left_variant'],c['right_variant'],c['budget'],c['endpoint'])==('amp80-200','cnn_diffusion',right,'full','full',20,endpoint))
                rows.append([tag,title,labels[right],f'{scale*c["mean_difference"]:+.3f}',f'[{scale*c["ci95"][0]:.3f}, {scale*c["ci95"][1]:.3f}]',str(c['n'])])
    tables+=tex_table('contrasts','Within-cohort, within-gene differences: DiffuPrimer minus comparator, with existing 95\\% paired bootstrap intervals. Count endpoints include all 75 genes; conditional endpoints use jointly defined genes. pp: percentage points.',
                     ['Cohort','Endpoint','Comparator','Difference','95\\% CI','$n$'],rows,'lllrrr')
    (OUT/'supplement_tables.tex').write_text(tables)
    # Keep the fully verified bibliography intact. Derive its compact print form.
    bibliography=OUT/'references.bib'
    if bibliography.exists():
        blocks=re.split(r'(?=^@)',bibliography.read_text(),flags=re.M)
        formatted=[]
        for block in blocks:
            def abbreviate(m):
                names=m.group(1).split(' and ')
                return '  author = {'+' and '.join(names[:3]+['others'] if len(names)>3 else names)+'},'
            block=re.sub(r'^  author = \{(.*)\},$',abbreviate,block,flags=re.M)
            block=block.replace('number = {suppl_1}','number = {Suppl. 1}')
            if re.search(r'^  doi = ',block,re.M) or block.startswith('@inproceedings'):
                block=re.sub(r'^  url = .*\n','',block,flags=re.M)
            if block.startswith('@article') or block.startswith('@inproceedings'):
                block=re.sub(r'^  doi = .*\n','',block,flags=re.M)
            block=block.replace('booktitle = {Advances in Neural Information Processing Systems}',
                                'booktitle = {{NeurIPS}}')
            formatted.append(block)
        (OUT/'submission_references.bib').write_text('% Derived print bibliography: first three authors + et al.; journal/conference metadata printed compactly. Full DOI and primary URLs remain in references.bib.\n'+''.join(formatted))
    manifest={'generated_at':datetime.now(timezone.utc).isoformat(),'primary_method':'cnn_diffusion','document_status':'AI-assisted research draft for author review; not asserted to comply with Bioinformatics LLM drafting policy',
      'source_summaries':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in source_paths],
      'result_values':records,'figures':[{'source':source,'sha256':sha(ROOT/source),'copy':f'figures/{name}.pdf'} for name,source in figures.items()],
      'main_figure_subsets':crop_records,
      'public_code_url':'','public_data_url':'','author_metadata':json.loads((ROOT/'paper/author_metadata.json').read_text())}
    if bibliography.exists():
        manifest['bibliography']={'verified_full_file':'references.bib','sha256':sha(bibliography),
           'print_file':'submission_references.bib','print_sha256':sha(OUT/'submission_references.bib'),
           'formatting':'Oxford SCIMED compact author lists; primary source URLs remain in full verified bibliography and verification records',
           'style_reference':'https://academic.oup.com/DocumentLibrary/Cybersecurity/Mini-style%20SciMed.pdf'}
    (OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'output':str(OUT),'macros':len(macros),'derived_records':len(records),'copied_figures':len(figures)}))


if __name__=='__main__':main()
