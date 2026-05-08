# Fig.3 Dedicated Workspace

这个文件夹是 Fig.3 数值审核/复现专用复制包。原始项目文件没有移动；
这里的内容是从主项目复制出来的快照，方便集中检查 Fig.3。

## 目录结构

- `mathematica/`
  - Wolfram/BHPToolkit 数值生成器。
  - 负责 Teukolsky radial、spheroidal/spherical angular functions、patch、densify、spherical SR。

- `src/`
  - Python 数值核心。
  - `physics_fig3.py` 负责读取 JSON、计算 plain sum / series reduction / Cesaro。

- `scripts/`
  - Python 审计和诊断脚本。
  - `audit_fig3_numeric_status.py` 是当前主审计入口。
  - `diag_fig3_w2_methods.py` 是历史 `Momega=2` 方法比较诊断。

- `data/`
  - Fig.3 相关 JSON 缓存。
  - 包含旧主缓存、patch 缓存、spherical SR 缓存和若干历史/备份缓存。

- `results/`
  - Fig.3 图片、metrics、audit JSON、历史输出。

- `docs/`
  - `fig3_solver_files_summary.md`：当前最重要的文件地图和数值状态说明。
  - `README_root_snapshot.md`：主项目 README 的复制快照。

## 当前推荐入口

在这个 `fig3_workspace/` 目录内运行：

```powershell
python reproduce.py fig3-audit
python reproduce.py fig3-methods
```

输出会写入本文件夹内的 `results/`。

## 当前主图

```text
results/Fig3_paper_methods_grid.png
```

对应 metrics：

```text
results/fig3_paper_methods_grid_metrics.json
results/fig3_numeric_audit.json
```

## 当前主数据

```text
data/fig3_data.json
data/fig3_patch_w2_a0_a099.json
data/fig3_sr_spherical.json
```

## 当前主求解器

```text
mathematica/generate_fig3_patch_w2_highspin.wls
mathematica/densify_fig3_patch_angles.wls
mathematica/generate_fig3_sr_spherical.wls
src/physics_fig3.py
main_fig3_paper_methods_grid.py
scripts/audit_fig3_numeric_status.py
```

## 审核提醒

当前 Fig.3 仍不是完全闭合的论文级数值结果：

- `Momega=0.6` 仍 fallback 到 legacy SR。
- `a=0.99, Momega=0.6` 目前 patch 覆盖到 `l=2..41`。
- `a=0.99, Momega=2` 目前 281 点角网格覆盖到 `l=2..29`，更高 `l` 仍保留旧短网格 radial 数据。
- `fig3_sr_spherical.json` 目前只含 `Momega=2` 的旧 spherical SR row，且 `l_mode_max=60`。

## 本轮修复入口

已经加入针对论文 Fig.3 的高频/高自旋修复流程：

```text
docs/FIG3_REPAIR_RUNBOOK.md
cluster/README_CLUSTER_FIG3.md
cluster/run_fig3_a099_w2_l30_70.slurm
cluster/run_fig3_a099_w06_l42_60.slurm
cluster/run_fig3_merge_sr_plot.slurm
```

修复目标：

- 补齐 `a=0.99, Momega=2, l=30..70` 的 281 点角网格 mode。
- 补齐 `a=0.99, Momega=0.6, l=42..60`。
- 重新生成 `fig3_sr_spherical.json`，默认 `FIG3_SR_LMAX=70`、`FIG3_SR_BASIS_PAD=24`。
- 重新审计、重画 `results/Fig3_paper_methods_grid.png`。
- 生成 `theta=120 deg` 的 `l_max` 收敛诊断。
