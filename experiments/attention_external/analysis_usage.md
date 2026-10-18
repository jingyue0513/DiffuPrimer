# 独立额外基因队列的分析

该队列单独使用 `data/processed/attention_external_targets.jsonl` 的 `external_test` 基因。不会载入原研究的 47 个结果来源，也不将原 75 个测试基因与额外 75 个基因合并。

```bash
.venv/bin/python experiments/attention/analyze.py \
  --targets data/processed/attention_external_targets.jsonl \
  --split external_test \
  --no-old-sources \
  --protocol experiments/attention_external/protocol.json \
  --evaluations experiments/attention_external/results/evaluations \
  --output experiments/attention_external/results/study_summary.json

.venv/bin/python experiments/attention/plot_results.py \
  --summary experiments/attention_external/results/study_summary.json \
  --output experiments/attention_external/results/figures \
  --reference-summary experiments/attention/results/study_summary.json
```

计算尚未全部完成时，在命令末尾加 `--allow-partial`。分析只接收完整基因队列、候选文件及已完成的完整参考评估，不把待筛查来源当成零结果。

默认完成分析要求 protocol 的 `expected_runs` 均完成。运行名、评价目录和主方法声明由实际实验协议提供；本说明不生成实验计划或结果。

输出为 `study_summary.json`、`.per_gene.jsonl`、`.groups.csv`、`.contrasts.csv`。指标、候选预算前缀与配对 bootstrap 使用相同计算函数；组内 CI 和方法比较以该队列基因为单位。每个组保留完整 N，包括没有返回候选的基因。

`05_external_generalization` 在两个队列存在已完成的共同方法时生成，展示分别计算的均值与区间；两个队列之间不配对不同基因。该队列来自相同物种与 RefSeq 来源，是额外留出基因评估，不是独立物种或新的实体检测数据。

可将独立结果追加到技术审阅网页的 3.1：

```bash
.venv/bin/python scripts/build_attention_research_page.py \
  --external-summary experiments/attention_external/results/study_summary.json
```

网页保持四节和八个实验问题。原队列与额外队列的表格、统计来源及数值验证标签分开；调试预览可以加 `--allow-partial`，会明确显示为技术审阅草稿。
