// Verify the paper website against recorded scientific data and actual browser rendering.
const fs=require('fs'),path=require('path');
const {chromium}=require(process.env.DIFFUPRIMER_PLAYWRIGHT_MODULE||'playwright');
(async()=>{
const root=path.resolve(__dirname,'..'),out=path.join(root,'reports/chinese_research');
const stats=JSON.parse(fs.readFileSync(path.join(root,'results/final_fixed_seed/study_summary.json'),'utf8'));
const values=JSON.parse(fs.readFileSync(path.join(root,'paper/reported_values.json'),'utf8')).values;
const b=await chromium.launch({headless:true,args:['--no-sandbox']});
const p=await b.newPage({viewport:{width:1440,height:1100},acceptDownloads:true});
const errors=[];p.on('pageerror',e=>errors.push(e.message));p.on('console',m=>{if(m.type()==='error')errors.push(m.text())});
const check=(x,message)=>{if(!x)throw new Error(message)};
await p.goto('http://127.0.0.1:8768/',{waitUntil:'load'});await p.waitForFunction(()=>window.researchReportReady===true);await p.evaluate(()=>document.fonts.ready);
check(await p.locator('main>section').count()===4,'Expected exactly four paper sections');
check(await p.locator('.contributions>li').count()===3,'Expected exactly three contributions');
check(await p.locator('.experiment').count()===8,'Expected eight experimental questions');
check(await p.locator('#method-figure').count()===1 && await p.locator('#method-figure img').evaluate(im=>im.complete&&im.naturalWidth>=3000),'English method figure must load at full resolution');
check(await p.locator('video').count()===0,'Video should not be present');
const text=await p.locator('main').innerText();check(!/264\s*组|369,000|运行进度|生成实验与完整参考库审计/.test(text),'Progress dashboard text remains');
const reported=await p.locator('[data-value]').evaluateAll(es=>es.map(e=>({key:e.dataset.value,digits:Number(e.dataset.digits||2),text:e.textContent})));
reported.forEach(r=>check(r.text===Number(values[r.key]).toFixed(r.digits),'Reported value differs: '+r.key));
const methods=['primer3','uniform_pair_sampler','rule_sampler','direct','diffusion'];
const group=(m,c='amp80-200',v='full',k=20)=>stats.summaries.find(r=>r.method===m&&r.context===c&&r.variant===v&&r.budget===k);
const plotted=async id=>JSON.parse(await p.locator('#'+id+' .plot').getAttribute('data-rows'));
for(const [id,key]of [['quality-count','unique_RNA_pass_pairs_without_genomic_flag'],['quality-fraction','RNA_no_genome_fraction']]){
  const rows=await plotted(id);check(rows.length===5,'Missing method in '+id);rows.forEach(r=>check(r.value===group(r.method)[key],'Chart differs: '+id+'/'+r.method));
}
for(const k of [1,2,5,10,20]){
 await p.locator('[data-budget="'+k+'"]').click();check(await p.locator('[data-budget="'+k+'"]').getAttribute('aria-pressed')==='true','Budget selection failed');
 const rows=await plotted('budget-count');rows.forEach(r=>check(r.value===group(r.method,'amp80-200','full',k).unique_RNA_pass_pairs,'Budget value differs'));
 check((await p.locator('#budget-unit').innerText()).includes('K = '+k),'Budget label differs');
}
for(const [id,contexts]of [['completion-count',['amp80-200_fix_forward','amp80-200_fix_reverse']],['redesign-count',['redesign_forward','redesign_reverse']]]){
 const rows=await plotted(id);rows.forEach((rs,i)=>rs.forEach(r=>check(r.value===group(r.method,contexts[i]).unique_RNA_pass_pairs,'Anchored chart differs')));
}
for(const [id,key]of [['steps-rna','unique_RNA_pass_pairs'],['steps-diversity','RNA_pass_amplicon_distance'],['steps-temperature','best_RNA_pass_temperature_penalty']]){
 const series=await plotted(id);series[0].values.forEach((v,i)=>check(v===group('diffusion','amp80-200',['steps1','steps4','steps8','full'][i])[key],'Denoising chart differs'));
}
const figureCount=await p.locator('figure.chart svg').count();check(figureCount===15,'Unexpected result plot count');
const links=await p.locator('.outline a,.mobile-nav a').evaluateAll(es=>es.map(e=>e.getAttribute('href').slice(1)));for(const id of links)check(await p.locator('[id="'+id+'"]').count()===1,'Broken section link '+id);
await p.locator('.outline a[href="#exp-redesign"]').click();await p.waitForFunction(()=>Math.abs(document.getElementById('exp-redesign').getBoundingClientRect().top-68)<5);check(await p.locator('.outline a[href="#exp-redesign"]').getAttribute('class')==='sub active','Active outline did not follow section');
const promise=p.waitForEvent('download');await p.locator('#quality-count .svg-download').click();const dl=await promise,file=path.join(out,'check_quality.svg');await dl.saveAs(file);const svg=fs.readFileSync(file,'utf8');check(svg.includes('7.36')&&svg.includes('<svg'),'Downloaded figure differs');fs.unlinkSync(file);
for(const name of ['manuscript.pdf','supplementary.pdf','study_summary.groups.csv','DiffuPrimer_method.pptx','DiffuPrimer_method.pdf','DiffuPrimer_method.svg','DiffuPrimer_method.png']){const response=await p.request.get('http://127.0.0.1:8768/'+name);check(response.status()===200,'Resource link failed '+name)}
await p.evaluate(()=>{window.scrollTo({top:0,behavior:'instant'});reading()});await p.screenshot({path:path.join(out,'preview_desktop.png')});
await p.locator('#exp-quality').scrollIntoViewIfNeeded();await p.screenshot({path:path.join(out,'preview_results.png')});
await p.setViewportSize({width:390,height:844});await p.evaluate(()=>{window.scrollTo({top:0,behavior:'instant'});reading()});check(await p.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),'Mobile overflow');await p.screenshot({path:path.join(out,'preview_mobile.png')});
await p.locator('#exp-quality').scrollIntoViewIfNeeded();check(await p.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),'Mobile chart overflow');await p.screenshot({path:path.join(out,'preview_mobile_results.png')});
const report={checked_at:new Date().toISOString(),presentation:'paper_website',sections:4,contributions:3,experiment_sections:8,method_figure_loaded:true,videos:0,result_plots:figureCount,budget_states:5,reported_value_instances:reported.length,exact_chart_values_verified:true,svg_download_verified:true,resource_links_verified:7,desktop_viewport:[1440,1100],mobile_viewport:[390,844],page_errors:errors,passed:errors.length===0};
fs.writeFileSync(path.join(out,'browser_check.json'),JSON.stringify(report,null,2)+'\n');await b.close();check(errors.length===0,'Browser errors: '+errors.join('; '));console.log(JSON.stringify(report,null,2));
})().catch(e=>{console.error(e);process.exit(1)});
