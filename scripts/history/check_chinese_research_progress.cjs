// Browser-level checks for the user-facing research report; no scientific data are simulated.
const fs = require('fs');
const path = require('path');
const { chromium } = require(process.env.DIFFUPRIMER_PLAYWRIGHT_MODULE || 'playwright');
(async () => {
  const root = path.resolve(__dirname, '..');
  const directory = path.join(root, 'reports/chinese_research');
  const browser = await chromium.launch({ headless: true, args: ['--no-sandbox'] });
  const page = await browser.newPage({ viewport: { width: 1440, height: 1100 }, acceptDownloads: true });
  const errors = [];
  page.on('pageerror', e => errors.push(e.message));
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text()); });
  await page.goto('file://' + path.join(directory, 'index.html'));
  await page.waitForFunction(() => window.researchReportReady === true);
  const check = (pass, message) => { if (!pass) throw new Error(message); };
  check(await page.locator('#training-table tbody tr').count() === 6, 'Missing training rows');
  check(await page.locator('#gene option').count() === 75, 'Missing test genes');
  for (const key of ['RNA_no_genome_fraction', 'unique_RNA_pass_pairs', 'two_RNA_pass_pairs', 'unique_RNA_pass_pairs_without_genomic_flag', 'two_RNA_pass_pairs_without_genomic_flag', 'RNA_pass_amplicon_distance', 'best_RNA_pass_temperature_penalty']) {
    await page.selectOption('#endpoint', key);
    for (const budget of ['1', '2', '5', '10', '20']) {
      await page.selectOption('#budget', budget);
      check(await page.locator('#main-table tbody tr').count() === 5, 'Missing main comparison method');
      check(!/NaN|undefined/.test(await page.locator('#main-table').innerText()), 'Invalid main result');
    }
  }
  await page.selectOption('#endpoint', 'unique_RNA_pass_pairs');
  await page.selectOption('#budget', '20');
  const statsPath=path.join(root,'results/final_fixed_seed',fs.existsSync(path.join(root,'results/final_fixed_seed/study_summary.json'))?'study_summary.json':'existing_quality.json');
  const stats=JSON.parse(fs.readFileSync(statsPath,'utf8'));
  const mainDiff=stats.summaries.find(r=>r.method==='diffusion'&&r.context==='amp80-200'&&r.variant==='full'&&r.budget===20);
  const mainContrast=stats.paired_comparisons.find(r=>r.context==='amp80-200'&&r.right==='primer3_full'&&r.endpoint==='unique_RNA_pass_pairs'&&r.budget===20);
  check((await page.locator('#main-table').innerText()).includes(mainDiff.unique_RNA_pass_pairs.toFixed(2)), 'Main diffusion number differs');
  check((await page.locator('#main-callout').innerText()).includes(mainContrast.mean_difference.toFixed(2)), 'Paired difference differs');
  for (const c of ['amp80-120','amp121-160','amp161-200','amp80-200_fix_forward','amp80-200_fix_reverse','redesign_forward','redesign_reverse']) {
    await page.selectOption('#condition', c);
    for (const metric of ['unique_RNA_pass_pairs','unique_RNA_pass_pairs_without_genomic_flag','RNA_pass_amplicon_distance','best_RNA_pass_temperature_penalty']) {
      await page.selectOption('#condition-metric', metric);
      check(await page.locator('#condition-chart svg').count() === 1, 'Condition chart absent');
    }
  }
  await page.selectOption('#condition', 'amp80-120');
  await page.selectOption('#condition-metric', 'unique_RNA_pass_pairs');
  // Exercise all genes and all methods, preserving every recorded no-design outcome.
  const geneIds = await page.locator('#gene option').evaluateAll(options => options.map(o => o.value));
  await page.evaluate(ids => {
    for (const id of ids) {
      document.getElementById('gene').value=id;
      for (const method of ['diffusion','direct','primer3','uniform_pair_sampler','rule_sampler']) {
        document.getElementById('candidate-method').value=method;
        document.getElementById('gene').dispatchEvent(new Event('change'));
        if(document.querySelectorAll('#gene-table tbody tr').length!==5)throw new Error('Gene statistics absent');
        if(document.querySelectorAll('#candidate-table tbody tr').length!==20)throw new Error('Recorded budget differs');
      }
    }
  }, geneIds);
  await page.fill('#gene-search', 'does-not-exist-xyz');
  check(await page.locator('#gene option').count() === 0, 'Empty search incorrect');
  await page.fill('#gene-search', '');
  await page.selectOption('#candidate-method','diffusion');
  await page.selectOption('#gene',geneIds[0]);
  await page.selectOption('#candidate-filter','nogdna');
  await page.selectOption('#candidate-filter','all');
  const figureButtons=page.locator('#figure-tabs button');
  check(await figureButtons.count() === 10, 'Missing figure');
  for(let i=0;i<10;i++) {
    await figureButtons.nth(i).click();
    await page.waitForFunction(()=>{let img=document.getElementById('figure-image');return img.complete&&img.naturalWidth>100});
  }
  await figureButtons.nth(0).click();
  for(const [button,name] of [['#export-main','comparison.csv'],['#export-candidates','candidates.csv'],['#export-summary','results.json']]) {
    const promise=page.waitForEvent('download');
    await page.locator(button).click();
    const download=await promise;
    const output=path.join(directory,'check_'+name);
    await download.saveAs(output);
    const contents=fs.readFileSync(output,'utf8');
    check(contents.length>100,'Empty download');
    if(name.endsWith('.json'))check(JSON.parse(contents).completion.completed_sampling_runs===264,'Exported statistics differ');
    fs.unlinkSync(output);
  }
  const links=await page.locator('a[href^="evidence/"]').evaluateAll(elements=>elements.map(a=>a.getAttribute('href')));
  for(const href of links)check(fs.existsSync(path.join(directory,href)),'Broken evidence link: '+href);
  await page.evaluate(()=>window.scrollTo(0,0));
  await page.screenshot({path:path.join(directory,'preview_desktop.png')});
  await page.setViewportSize({width:390,height:844});
  await page.evaluate(()=>window.scrollTo(0,0));
  check(await page.evaluate(()=>document.documentElement.scrollWidth<=window.innerWidth+1),'Mobile page overflows');
  await page.screenshot({path:path.join(directory,'preview_mobile.png')});
  const report={checked_at:new Date().toISOString(),browser:'Chromium',page_errors:errors,main_comparison_states:35,condition_states:28,gene_method_states:375,figures:10,downloads_checked:3,evidence_links_checked:links.length,desktop_viewport:[1440,1100],mobile_viewport:[390,844],passed:errors.length===0};
  fs.writeFileSync(path.join(directory,'browser_check.json'),JSON.stringify(report,null,2)+'\n');
  await browser.close();
  check(errors.length===0,'Browser errors: '+errors.join('; '));
  console.log(JSON.stringify(report,null,2));
})().catch(e=>{console.error(e);process.exit(1)});
