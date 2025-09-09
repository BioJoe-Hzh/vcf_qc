#!/usr/bin/env python3
"""
演示 VCF QC Site Level 新功能的测试脚本
展示详细的过滤统计报告和新的阈值控制
"""

import pandas as pd
import numpy as np

# 模拟创建测试数据
def create_test_data():
    """创建模拟的site metrics数据用于演示"""
    np.random.seed(42)
    n_sites = 1000
    
    data = {
        'Chrom': ['chr1'] * n_sites,
        'Pos': range(1, n_sites + 1),
        'QUAL': np.random.exponential(25, n_sites),  # 指数分布，一些低质量位点
        'MeanDepth': np.random.gamma(20, 2, n_sites),  # 深度分布
        'QD': np.random.gamma(8, 2, n_sites),  # QD分布，一些异常值
        'MAC': np.random.poisson(5, n_sites),  # MAC分布，包含单例变异
        'MAF': np.random.beta(1, 20, n_sites),  # MAF分布，大多数低频
        'MissingRate': np.random.beta(1, 50, n_sites),  # 大多数位点missing rate很低
    }
    
    # 添加一些异常值来展示过滤效果
    data['QUAL'][:50] = np.random.uniform(1, 15, 50)  # 一些低质量位点
    data['QD'][:30] = np.random.uniform(0.5, 1.5, 30)  # 一些低QD位点
    data['QD'][-20:] = np.random.uniform(40, 60, 20)  # 一些异常高QD位点
    data['MissingRate'][-100:] = np.random.uniform(0.2, 0.8, 100)  # 一些高missing rate位点
    data['MAC'][:200] = 1  # 一些单例变异
    
    return pd.DataFrame(data)

if __name__ == "__main__":
    print("🧪 VCF QC Site Level 新功能演示")
    print("=" * 50)
    
    # 创建测试数据
    df = create_test_data()
    print(f"📊 创建了包含 {len(df):,} 个位点的测试数据")
    
    # 显示数据分布概况
    print("\n📈 数据分布概况:")
    for col in ['QUAL', 'QD', 'MAC', 'MAF', 'MissingRate']:
        values = df[col]
        print(f"   {col}: min={values.min():.3f}, max={values.max():.3f}, mean={values.mean():.3f}")
    
    print("\n" + "=" * 50)
    print("🔍 以下是使用新功能时的预期输出示例:")
    print("=" * 50)
    
    print("""
命令: vcf_qc site --vcf test.vcf.gz --out results \\
        --min-qual 20 --min-qd 2.0 --max-qd 35.0 \\
        --max-missing 0.1 --min-mac 2

预期输出:
🔍 Site-level filtering report:
   Initial sites: 1,000
   QUAL >= 20: removed 45 sites, 955 remaining
   QD >= 2.0: removed 28 sites, 927 remaining  
   QD <= 35.0: removed 18 sites, 909 remaining
   Missing rate <= 0.1: removed 95 sites, 814 remaining
   MAC >= 2: removed 180 sites, 634 remaining

📊 Final filtering summary:
   Total sites removed: 366 (36.6%)
   Sites retained: 634 (63.4%)
    """)
    
    print("=" * 50)
    print("✨ 新功能亮点:")
    print("  1. 📊 详细的逐步过滤统计")
    print("  2. 🎯 QD最小值和最大值控制")
    print("  3. 📈 每一步的保留率信息")
    print("  4. 🔍 清晰的过滤效果展示")
    print("=" * 50)
