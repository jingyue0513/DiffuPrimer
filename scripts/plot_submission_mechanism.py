#!/usr/bin/env python3
"""Render Supplementary Figure S5 from frozen mechanism summaries."""
import hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'paper/bioinformatics_submission'
def main():
    source=OUT/'mechanism_evidence.json';data=json.loads(source.read_text())
    methods=['direct_validation_selected','diffusion1_validation_selected','diffusion_4','frozen_4','diffusion_16','frozen_16','direct_equal_cpu_uncapped']
    labels=['Direct\n(validation-selected)','One-step denoising\n(validation-selected)','Diffusion, 4 steps','Frozen scores, 4 steps','Diffusion, 16 steps','Frozen scores, 16 steps','Direct, matched CPU']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.2,'pdf.fonttype':42,'axes.linewidth':.6})
    fig,axes=plt.subplots(1,2,figsize=(178/25.4,108/25.4),sharey=True)
    used=[]
    for ax,cohort,title in zip(axes,['test','external_test'],['A  Original cohort (75 genes)','B  Additional cohort (75 genes)']):
        gs=[next(g for g in data['groups'] if g['cohort']==cohort and g['method']==m) for m in methods];used.extend(gs)
        v=np.array([g['retained_pairs'] for g in gs]);ci=np.array([g['retained_pairs_ci95'] for g in gs]);y=np.arange(len(methods))
        colors=['#9174ad' if 'validation_selected' in m else '#dca04d' if 'frozen' in m else '#b9697b' if m=='direct_equal_cpu_uncapped' else '#5185a0' for m in methods]
        ax.barh(y,v,color=colors,height=.63,xerr=np.maximum(0,[v-ci[:,0],ci[:,1]-v]),error_kw={'elinewidth':.7,'capsize':2,'capthick':.7})
        for yy,value,hi in zip(y,v,ci[:,1]):ax.text(hi+.16,yy,f'{value:.2f}',va='center',fontsize=7)
        ax.set_xlim(0,14);ax.set_xticks([0,4,8,12]);ax.set_title(title,fontsize=8,loc='left',pad=9)
        ax.set_xlabel('Retained pairs per gene',fontsize=7.5);ax.spines[['top','right']].set_visible(False)
        ax.set_yticks(y,labels);ax.tick_params(axis='y',length=0,pad=5);ax.tick_params(axis='x',width=.6,length=3)
    axes[0].invert_yaxis();axes[1].tick_params(axis='y',labelleft=False)
    fig.subplots_adjust(left=.28,right=.99,top=.91,bottom=.12,wspace=.16)
    fig.savefig(OUT/'figures/mechanism_controls.pdf');fig.savefig(OUT/'figures/mechanism_controls.png',dpi=220);plt.close(fig)
    (OUT/'mechanism_figure_provenance.json').write_text(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'groups':used,'methods':methods,'display_labels':labels},indent=2)+'\n')
if __name__=='__main__':main()
