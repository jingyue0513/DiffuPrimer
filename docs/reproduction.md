# 复现与执行

正式版本采用 seed 20261013、现有 30 轮检查点及验证集选择的温度 4。最终分析包含改变长度条件后的重设计与反向锚点修正。使用已完成候选结果重建论文和中文网页：

```bash
PYTHONPATH=src .venv/bin/python scripts/analyze_fixed_seed_study.py
PYTHONPATH=src .venv/bin/python scripts/plot_fixed_seed_study.py
PYTHONPATH=src .venv/bin/python scripts/finalize_paper.py
PYTHONPATH=src .venv/bin/python scripts/verify_final_study.py
PYTHONPATH=src .venv/bin/python scripts/build_chinese_research_page.py
```

补充实验驱动：`experiments/final_followup/run.py` 与 `experiments/final_followup/run_anchor_prefilter.py`；已完成任务会复用落盘记录。正式分析替换四组受反向锚点修正影响的神经方法输出；原结果仍保留。协议见 `experiments/final_followup/protocol.json` 和 `anchor_prefilter_protocol.json`。完整比较表与每基因记录见 `results/final_fixed_seed/study_summary.*`。历史多种子实验无需重新运行。

以下保留原研究各阶段的执行命令。

工作目录为 `/share/project/yenaisheng/DiffuPrimer`。当前运行使用本地 `.venv`、CPU PyTorch 和每个计算进程 8 个线程。`requirements.lock` 保存实际安装版本；MMseqs2 与 BLAST+ 二进制由官方来源获取，版本与文件校验值另行记录。GPU 上已有其他任务，当前研究使用 CPU。

参考数据固定为 GCF_000001635.27 / GCF_000001635.27-RS_2024_02：

```bash
PYTHONPATH=src .venv/bin/python scripts/download_reference.py --genome
PYTHONPATH=src .venv/bin/python scripts/verify_reference.py --genome
PYTHONPATH=src .venv/bin/python scripts/build_specificity_reference.py --genome
```

本次研究已冻结 500 个基因抽样及训练/验证/测试划分。现有来源访问被拒绝后，新采集不会自动重试；缓存可继续使用。`scripts/acquire_pilot.py` 在服务拒绝或限流时终止新请求。

教师候选只由训练和验证基因的固定模板生成，所有结果明确为计算候选：

```bash
PYTHONPATH=src .venv/bin/python scripts/generate_weak_supervision.py
PYTHONPATH=src .venv/bin/python scripts/prepare_data.py --freeze-partial --include-teacher
PYTHONPATH=src .venv/bin/python scripts/freeze_protocol.py
PYTHONPATH=src .venv/bin/python scripts/run_study.py --extended
```

`run_study.py` 运行三种子的直接 CNN 和扩散 CNN、匹配候选预算的规则与 Primer3 对照、验证集、完整测试集、RNA/基因组筛查、汇总与论文构建。扩展运行包含长度条件和单侧固定补全。每条命令的输出保存在 `results/study_logs/`。

校准检查只运行验证集，输出目录带 `calibration`。这些结果不会进入论文正式比较表。完整测试集没有按可设计性删去目标。已冻结代码或数据变化后，必须记录版本更正并完整重算受到影响的结果，不能静默混用旧结果。

```bash
PYTHONPATH=src .venv/bin/python scripts/build_paper.py --aggregate results/primary_aggregate.json
```

`build_paper.py` 现在调用正式论文构建入口 `finalize_paper.py`，同时生成正文和补充材料。论文研究范围与用户最终确认一致，为公开数据与系统计算实验。
