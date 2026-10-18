// Check the separate technical review locally, including scientific values/offline assets.
const fs = require('fs'), path = require('path'), crypto = require('crypto'), {pathToFileURL} = require('url');
const {chromium} = require(process.env.DIFFUPRIMER_PLAYWRIGHT_MODULE || 'playwright');
(async () => {
  const root = path.resolve(__dirname, '..'), out = path.join(root, 'reports/cnn_publication');
  const summary = JSON.parse(fs.readFileSync(path.join(out, 'study_summary.json'), 'utf8'));
  const groups = new Map(summary.summaries.map(g => [[g.context,g.method,g.variant,g.budget].join('|'),g]));
  const externalPath=path.join(out,'external_study_summary.json');
  const external=fs.existsSync(externalPath)?JSON.parse(fs.readFileSync(externalPath,'utf8')):null;
  const externalGroups=new Map(external?external.summaries.map(g=>[[g.context,g.method,g.variant,g.budget].join('|'),g]):[]);
  const assert = (value,message) => {if (!value) throw new Error(message);};
  const browser = await chromium.launch({headless:true,args:['--no-sandbox']});
  const page = await browser.newPage({viewport:{width:1440,height:1100}});
  const errors = [], remoteRequests = [];
  page.on('pageerror', e => errors.push(e.message));
  page.on('request', r => {if (/^https?:/.test(r.url())) remoteRequests.push(r.url());});
  await page.goto(pathToFileURL(path.join(out,'index.html')).href,{waitUntil:'load'});
  await page.waitForFunction(() => window.researchReportReady === true);
  await page.evaluate(() => document.fonts.ready);
  assert(await page.locator('main>section').count() === 4, 'Need four primary sections');
  assert(await page.locator('.contributions>li').count() === 3, 'Need three separate contributions');
  assert(await page.locator('.experiment').count() === 8, 'Need eight experimental questions');
  assert(await page.locator('video').count() === 0, 'Video should not be present');
  const authorMetadata=JSON.parse(fs.readFileSync(path.join(root,'paper/author_metadata.json'),'utf8'));
  const authorText=await page.locator('.author-line').innerText();
  for (const author of authorMetadata.authors) {
    assert(authorText.includes(author.name) && authorText.includes(author.institution), 'Author metadata missing');
    if (author.corresponding && author.email) assert(await page.locator(`.author-line a[href="mailto:${author.email}"]`).count()===1, 'Corresponding email missing');
  }
  for (const figure of await page.locator('.result-figure img').all()) {
    await figure.scrollIntoViewIfNeeded();
    await figure.evaluate(im => im.decode());
    assert(await figure.evaluate(im => im.complete && im.naturalWidth >= 2000), 'Publication figure failed to load');
  }
  assert(await page.locator('.result-figure').count() === 5, 'Need five integrated figures');
  assert(await page.locator('#method-figure img').evaluate(im => im.complete && im.naturalWidth >= 2000), 'Biological CNN figure must load at publication resolution');
  const checkCells = async () => {
    const cells = await page.locator('td[data-endpoint]').evaluateAll(es => es.map(e => ({...e.dataset})));
    for (const cell of cells) {
      const lookup=cell.cohort==='external'?externalGroups:groups;
      const g = lookup.get([cell.context,cell.method,cell.variant,Number(cell.budget)].join('|'));
      assert(g && JSON.parse(cell.value) === g[cell.endpoint], 'Website data diverged: ' + JSON.stringify(cell));
    }
    return cells.length;
  };
  const initialCells = await checkCells();
  for (const k of [1,2,5,10,20]) {
    await page.locator(`.budget-buttons button[data-budget="${k}"]`).click();
    assert(await page.locator(`.budget-buttons button[data-budget="${k}"]`).getAttribute('aria-pressed') === 'true', 'Budget selection failed');
    assert(await page.locator('#budget-results td[data-budget]').evaluateAll((es,k) => es.every(e => Number(e.dataset.budget) === k), k), 'Budget cell attributes stale');
    await checkCells();
  }
  const text = await page.locator('main').innerText();
  assert(!/264\s*组|运行进度|完成率|训练进度/.test(text), 'Progress dashboard language remains');
  if (summary.status !== 'complete') assert(await page.locator('.draft-label').count() === 1, 'Working preview must be explicit');
  const links = await page.locator('.outline a,.mobile-nav a').evaluateAll(es => es.map(e => e.getAttribute('href').slice(1)));
  for (const id of links) assert(await page.locator(`[id="${id}"]`).count() === 1, 'Broken section target');
  const resources = await page.locator('a[href]').evaluateAll(es => es.map(e => e.getAttribute('href')).filter(h => !h.startsWith('#') && !/^(https?:|mailto:)/.test(h)));
  for (const resource of resources) assert(fs.existsSync(path.join(out,resource)), 'Missing download: ' + resource);
  await page.evaluate(() => window.scrollTo({top:0,behavior:'instant'}));
  await page.screenshot({path:path.join(out,'preview_desktop.png')});
  await page.locator('#exp-quality').scrollIntoViewIfNeeded();
  await page.screenshot({path:path.join(out,'preview_results.png')});
  await page.setViewportSize({width:390,height:844});
  await page.evaluate(() => window.scrollTo({top:0,behavior:'instant'}));
  assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1), 'Mobile horizontal overflow');
  await page.screenshot({path:path.join(out,'preview_mobile.png')});
  await page.locator('#exp-quality').scrollIntoViewIfNeeded();
  assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1), 'Result table causes mobile overflow');
  await page.screenshot({path:path.join(out,'preview_mobile_results.png')});
  assert(errors.length === 0, 'Browser script errors: ' + errors.join('; '));
  assert(remoteRequests.length === 0, 'Offline page requests remote resources');
  const report = {at:new Date().toISOString(),status:summary.status,passed:true,
    page_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(out,'index.html'))).digest('hex'),
    summary_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(out,'study_summary.json'))).digest('hex'),
    sections:4,contributions:3,experimental_questions:8,
    numeric_cells_verified:initialCells,budget_states_verified:5,method_figure_loaded:true,offline_resources:true,
    author_metadata_verified:true,
    separate_external_cohort:external?{total_genes:external.total_test_genes,split:external.evaluation_split,old_sources_included:external.old_sources_included}:null,
    desktop_viewport:[1440,1100],mobile_viewport:[390,844],resource_links_verified:resources.length,page_errors:errors,remote_requests:remoteRequests};
  fs.writeFileSync(path.join(out,'browser_check.json'),JSON.stringify(report,null,2)+'\n');
  await browser.close();
  console.log(JSON.stringify(report,null,2));
})().catch(error => {console.error(error);process.exit(1);});
