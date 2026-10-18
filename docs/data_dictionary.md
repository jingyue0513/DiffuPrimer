# 数据与证据规则

一条记录代表一对引物和一个明确版本的目标模板。`gene_id` 是返回并核对过的 NCBI GeneID，`transcript_accession_version` 是当前参考中的带版本 accession。`source_transcript_accession`、`source_validation_accession_version` 和 `source_template_sequence` 保留原数据库信息，不能直接替换为当前参考。

`forward_start/end` 与 `reverse_start/end` 均为模板正链上的 **0-based、半开区间**；引物之间不重叠。`reverse_primer` 是实际合成方向的 5′→3′ 反向互补序列。`amplicon_length = reverse_end - forward_start`，包含两条引物。

`mapping_status` 保存唯一精确定位、源模板歧义、当前参考无法定位或基因解析不符等结果。位置重新映射不等于在新参考或新样本上重新完成实验验证。每次原始 HTTP 响应保存 URL、查询参数、时间和 SHA-256。

`evidence_type` 有四类：`experimental_success`、`experimental_failure`、`computational_candidate`、`unverified`。`experimental_record_available` 只表示存在实验材料。没有成功判定的曲线、BLAST 页面、绿色链接或用户反馈表单不能赋成功标签；未知、未测试、实验失败和设计失败彼此区分。公开论文里的明确结果必须同时核对 PrimerBank ID、两条序列和产物长度才能关联。

Primer3 生成的数据属于 `computational_candidate`，其 `experimental_result` 为 `not_tested`。这类数据可训练规则模仿模型，不能统计为湿实验成功。公开来源不足导致的部分数据训练需保存缺失/错误清单，并明确标注弱监督研究阶段。

`geometry_pass`、`thermodynamic_pass`、`computational_feasible` 是计算检查。`gene_specificity_screen_pass` 仅表示通过当前全长无 gap、≤2 错配、3′端 5 nt 完全匹配和 36–1000 bp 扩增范围的局部 BLAST 规则。未知基因注释防止筛查判定为通过。`genomic_DNA_product_detected` 单独报告潜在 gDNA 扩增。两者均不能替代独立实验或完整 Primer-BLAST 复核。

每个基因以 K=20 个候选预算评估；无效输出、重复输出和 Primer3 返回不足均不补抽。唯一引物序列对、产物位置差异、产物交并比距离需分别报告，避免把 1 nt 位置变化夸大为大量独立方案。统计单位是配对的基因，候选与种子不作为额外生物重复。
