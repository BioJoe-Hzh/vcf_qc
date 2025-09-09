# VCF QC Site Level 改进功能使用指南

## 新增功能概览

### 1. 全局Smart Cutoff控制
- `--no-smart-cutoff`: 禁用所有指标的99.5%智能截断
- 显示数据的完整分布范围

### 2. 特定指标Smart Cutoff控制
- `--no-cutoff-metrics`: 只对指定指标禁用smart cutoff
- 支持的指标：`qual`, `depth`, `qd`, `mac`, `maf`, `missing`
- 可以指定多个指标，用空格分隔

### 3. 扩展的过滤阈值控制
- `--min-qual`: QUAL最小阈值
- `--max-missing`: Missing rate最大阈值  
- `--min-depth`: 深度最小阈值
- `--max-depth`: 深度最大阈值 (新增)
- `--min-mac`: MAC最小阈值
- `--min-qd`: QD最小阈值

### 4. MAC/MAF特定阈值过滤
- `--mac-plot-min`: 设置MAC绘图的最小值阈值
- `--maf-plot-min`: 设置MAF绘图的最小值阈值

### 5. 详细过滤统计报告 (新增)
- 自动显示每一步过滤的效果
- 显示每个阈值过滤掉的位点数量
- 显示最终保留的位点数量和比例

## 使用示例

### 基础使用（保持原有行为）
```bash
vcf_qc site --vcf input.vcf.gz --out site_results
```

### 质量控制过滤示例
```bash
# 基础质控过滤
vcf_qc site --vcf input.vcf.gz --out site_results \
    --min-qual 30 \
    --min-qd 2.0 \
    --max-depth 1000 \
    --max-missing 0.1 \
    --min-mac 2

# 输出示例：
# 🔍 Site-level filtering report:
#    Initial sites: 100,000
#    QUAL >= 30: removed 5,234 sites, 94,766 remaining
#    Missing rate <= 0.1: removed 3,412 sites, 91,354 remaining
#    Mean depth <= 1000: removed 1,234 sites, 90,120 remaining  
#    MAC >= 2: removed 12,345 sites, 77,775 remaining
#    QD >= 2.0: removed 2,145 sites, 75,630 remaining
#
# 📊 Final filtering summary:
#    Total sites removed: 24,370 (24.4%)
#    Sites retained: 75,630 (75.6%)
```

### 组合使用：过滤 + 特定显示控制
```bash
# 严格质控 + MAF完整分布显示
vcf_qc site --vcf input.vcf.gz --out site_results \
    --min-qual 30 --min-qd 2.0 --max-missing 0.05 \
    --no-cutoff-metrics maf missing \
    --maf-plot-min 0.01
```

### 探索性分析（显示完整分布）
```bash
vcf_qc site --vcf input.vcf.gz --out site_results \
    --no-cutoff-metrics maf missing qd
```

## 参数详解

### 过滤阈值参数
| 参数 | 类型 | 说明 | 推荐值 |
|------|------|------|--------|
| `--min-qual` | float | QUAL最小值 | 20-50 |
| `--max-missing` | float | Missing rate最大值 | 0.05-0.2 |
| `--min-depth` | int | 平均深度最小值 | 8-20 |
| `--max-depth` | int | 平均深度最大值 | 1000-2000 |
| `--min-mac` | int | MAC最小值 | 1-5 |
| `--min-qd` | float | QD最小值 | 2.0-5.0 |

### 显示控制参数
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--no-smart-cutoff` | flag | False | 禁用所有指标的smart cutoff |
| `--no-cutoff-metrics` | list | [] | 只对指定指标禁用smart cutoff |
| `--mac-plot-min` | int | 0 | MAC绘图最小阈值 |
| `--maf-plot-min` | float | 0.0 | MAF绘图最小阈值 |

## 过滤统计报告解读

### 报告格式示例
```
🔍 Site-level filtering report:
   Initial sites: 100,000
   QUAL >= 30: removed 5,234 sites, 94,766 remaining
   QD >= 2.0: removed 2,145 sites, 92,621 remaining
   Missing rate <= 0.1: removed 3,412 sites, 89,120 remaining
   MAC >= 2: removed 12,345 sites, 76,775 remaining

📊 Final filtering summary:
   Total sites removed: 23,225 (23.2%)
   Sites retained: 76,775 (76.8%)
```

### 关键指标解释
- **Initial sites**: 原始位点总数
- **removed X sites**: 当前过滤步骤移除的位点数
- **X remaining**: 当前步骤后剩余的位点数
- **Total sites removed**: 所有过滤步骤累计移除的位点数
- **Sites retained**: 最终保留的位点数和比例

## 应用场景

### 1. 标准质控流程
```bash
vcf_qc site --vcf input.vcf.gz --out results \
    --min-qual 30 \
    --min-qd 2.0 \
    --max-missing 0.1 \
    --min-mac 2
```
- 适用于：常规变异调用质控

### 2. 严格质控 + 完整分析
```bash
vcf_qc site --vcf input.vcf.gz --out results \
    --min-qual 50 --min-qd 5.0 --max-missing 0.05 \
    --no-cutoff-metrics maf missing \
    --maf-plot-min 0.01
```
- 适用于：高质量数据集分析

### 3. 探索性质控评估
```bash
vcf_qc site --vcf input.vcf.gz --out results \
    --no-cutoff-metrics qd maf missing
```
- 适用于：数据质量初步评估

### 4. 发表级质控
```bash
vcf_qc site --vcf input.vcf.gz --out results \
    --min-qual 30 --min-qd 2.0 --max-depth 1000 \
    --max-missing 0.05 --min-mac 3 \
    --mac-plot-min 1 --maf-plot-min 0.01
```
- 适用于：发表文章的严格质控

## 技术细节

### 过滤顺序
1. QUAL过滤
2. Missing rate过滤  
3. 深度过滤
4. MAC过滤
5. QD最小值过滤
6. QD最大值过滤

### 输出文件
- `site_metrics_raw.tsv`: 原始指标数据
- `site_metrics_filtered.tsv`: 过滤后的指标数据
- 各种图表文件（PNG格式）

### 性能特点
- 过滤是逐步进行的，可以清楚看到每一步的效果
- 统计报告不影响处理速度
- 支持大规模数据集处理
