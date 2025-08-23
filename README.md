## vcf_qc – VCF 质量控制与可视化工具包

面向中小规模到中等规模 VCF 的轻量级 QC 指标与绘图。提供 sample / site / genotype 三层级指标与分布、联合散点、直方图、等位基因组成等，并支持多进程绘图加速。

### 主要特性
* 轻量级行流式 VCF 读取（不依赖 cyvcf2，纯 Python 解析基本字段）。
* Sample 级：缺失率、平均深度、杂合比 (Het/HomAlt)、自动/手动聚类 (可选 scikit-learn)。
* Site 级：QUAL、MeanDepth、QD、MAC、MAF(定义为第二多等位基因频率)、MissingRate、重构多等位基因位点的 AlleleCount 组成饼图。
* Genotype 级：GQ、Depth、VAF (多等位基因感知)、ALT_Read_Count，杂合子专用散点。
* 智能尾部截断 (smart cutoff)：迭代裁剪极端长尾，突出主分布。
* 多进程并行绘图：--threads N。
* 可选安装 `vcf_qc[cluster]` 激活聚类（依赖 scikit-learn）。

### 目录结构（简化）
```
vcf_qc/
  cli.py                 # 命令行入口逻辑（sample/site/genotype）
  io/                    # SimpleVCFReader
  metrics/               # sample_metrics / site_metrics 等
  plot/                  # 基础绘图与各级别 plot_* 函数
pyproject.toml           # 包配置 (PEP 621)
README.md
LICENSE
```

### 安装
开发模式（建议先在虚拟环境中）：
```bash
pip install -e .
```
带聚类功能：
```bash
pip install -e .[cluster]
```

未来发布到 PyPI 后可：
```bash
pip install vcf_qc
```

### 命令行使用
```
vcf_qc sample   --vcf input.vcf.gz --out sample_qc   --threads 4
vcf_qc site     --vcf input.vcf.gz --out site_qc     --threads 4 --min-qual 30 --min-mac 2
vcf_qc genotype --vcf input.vcf.gz --out genotype_qc --threads 4
```
核心通用参数：
* `--vcf`：输入 (bgzip/gzip/plain) VCF。
* `--out`：输出目录（自动创建）。
* `--max-site`：限制解析的位点数（调试/采样）。
* `--threads`：并行绘图进程数。

Site 子命令附加过滤：`--min-qual --max-missing --min-depth --min-mac --focus-qual-qd-pct`。

### 主要输出示例
Sample 级：
* sample_missing_rate.png / sample_mean_depth.png / sample_het_ratio.png
* depth_vs_missing_rate.png / het_ratio_vs_depth.png
* gq_boxplot_vs_depth.png

Site 级：
* site_qual_distribution.png / site_mean_depth_distribution.png
* site_qual_vs_mean_depth.png / site_qd_distribution.png
* site_mac_vs_mean_depth.png / site_maf_vs_mean_depth.png
* site_mac_distribution.png / site_maf_distribution.png
* site_missing_rate_distribution.png
* site_allele_count_pie.png

Genotype 级：
* genotype_gq_distribution.png / genotype_depth_distribution.png
* genotype_gq_vs_depth.png
* genotype_vaf_vs_depth_het.png / genotype_alt_count_vs_depth_het.png

### Python 内部调用示例
```python
from vcf_qc.io import SimpleVCFReader
from vcf_qc.metrics.sample_metrics import compute_sample_metrics
from vcf_qc.plot.sample_plots import plot_missing_rate_per_sample

reader = SimpleVCFReader("input.vcf.gz", max_records=50000)
df = compute_sample_metrics(reader)
plot_missing_rate_per_sample(dict(zip(df.Sample, df.MissingRate)), output_path="missing.png")
```

### 可选：在 README 中显示版本
```bash
python -c "import vcf_qc; print(vcf_qc.__version__)"
```

### 发布到 GitHub & 版本标签
1. 初始化 Git（若未初始化）：
    ```bash
    git init
    git add .
    git commit -m "feat: initial public release v0.1.0"
    ```
2. 创建 GitHub 公开仓库（Web 上 New Repository，名称建议：`vcf_qc`）。
3. 关联远程并推送：
    ```bash
    git remote add origin git@github.com:<USERNAME>/vcf_qc.git
    git branch -M main
    git push -u origin main
    ```
4. 打标签：
    ```bash
    git tag -a v0.1.0 -m "vcf_qc 0.1.0 initial release"
    git push origin v0.1.0
    ```

### （可选）发布到 PyPI
```bash
python -m pip install --upgrade build twine
python -m build
twine upload dist/*
```

### 许可证
MIT License (见 LICENSE)。

### 贡献指南（简略）
* Fork & PR。
* 代码格式：建议 ruff / black（可在未来添加）。
* 新增绘图请保持函数签名与现有 `plot_*` 风格一致，提供 docstring。

---
欢迎 Issue / PR 提交改进建议。

### 作者 / 联系方式
Zihao Huang  
Department of Genetics, University of Cambridge  
Email: zh384@cam.ac.uk  
