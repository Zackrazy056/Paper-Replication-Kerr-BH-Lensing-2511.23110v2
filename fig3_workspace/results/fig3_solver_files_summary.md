# Fig.3 数值求解器相关文件最新汇总

更新时间：2026-05-08

## 一、核心数值链路

Fig.3 当前链路分成三层：

1. Wolfram/Mathematica 生成 Teukolsky 模态数据：
   `Re_C`、`B_ratio`、`S_0`、`S_theta`、`S_pi_minus_theta`。

2. Python `Fig3Engine` 读取缓存 JSON，计算散射截面：
   `|f|^2 + |g|^2`，并实现 series reduction / Cesaro。

3. Python 画图与审计：
   生成 2x2 三方法图，并输出覆盖率、截断、raw max 等诊断。

## 二、Wolfram 求解器文件

### `mathematica/generate_fig3_patch_w2_highspin.wls`

当前 Fig.3 paper-grid patch 主求解器。

默认目标：

- `a = 0, 0.99`
- `Momega = 0.6, 2.0`
- `l = 2..90`
- `theta = 40..180 deg`，步长 `0.5 deg`，共 281 点
- 输出：`data/fig3_patch_w2_a0_a099.json`

当前重要修正：

- `a=0` 时直接用 `SpinWeightedSphericalHarmonicY`，避免 spheroidal 展开卡死。
- 修复 `ComputeHarmonicsWithPrec` 里 Mathematica `Module` 局部变量初始化问题。
- radial 尝试顺序改为先用带 `AccuracyGoal/PrecisionGoal` 的路径，再 fallback 到 no-goal。
- `FIG3_PATCH_CHECK_BRATIO` 默认关闭，避免高精度二次 bratio 校验阻塞整批。
- 旧的 71 点短角度行不会再被当作完整完成行。

常用命令：

```powershell
$env:FIG3_PATCH_BATCH = "20"
wolframscript -file mathematica/generate_fig3_patch_w2_highspin.wls
```

可用环境变量：

- `FIG3_PATCH_A_LIST`
- `FIG3_PATCH_OMEGA_LIST`
- `FIG3_PATCH_LMIN`
- `FIG3_PATCH_LMAX`
- `FIG3_PATCH_BATCH`
- `FIG3_PATCH_OUT`
- `FIG3_PATCH_CHECK_BRATIO`

### `mathematica/densify_fig3_patch_angles.wls`

角度网格补密工具。

用途：不重算 `B_ratio` / radial，只对已有短网格行重算
`S_0`、`S_theta`、`S_pi_minus_theta` 到 281 点。

默认目标：

- 输入/输出：`data/fig3_patch_w2_a0_a099.json`
- 默认 densify：`a=0.99, Momega=2`

常用命令：

```powershell
$env:FIG3_DENSIFY_A_LIST = "0.99"
$env:FIG3_DENSIFY_OMEGA_LIST = "2"
wolframscript -file mathematica/densify_fig3_patch_angles.wls
```

注意：当前 `a=0.99, Momega=2` 只成功 densify 到 `l=2..29`；
`l>=30` 的角函数本地计算仍很慢或失败，脚本会保留原行并标记
`densify_failed`，不会再破坏已有 radial 数据。

### `mathematica/generate_fig3_sr_spherical.wls`

Appendix-C aligned spherical-basis series-reduction 缓存生成器。

用途：

- 从 spheroidal modal `f(theta)` 投影到 spin-weighted spherical basis。
- 对 spherical-basis 系数做 series reduction。
- 输出预计算的 `f_sr(theta)` 到 `data/fig3_sr_spherical.json`。

当前状态：

- 已加入 `Momega=0.6` 的 combo target。
- 目前缓存仍只有两行：
  `a=0, Momega=2` 与 `a=0.99, Momega=2`。
- 当前缓存 `l_mode_max=60`，不是完整 `l=88/90`。
- 本地重建全 combo 时曾超过 30 分钟未落盘，仍是瓶颈。

### 旧版/历史生成器

`mathematica/generate_fig3_data.wls`

- 旧主数据生成器。
- 默认 `Momega=0.5,2.0`，不是新图要求的 `0.6,2.0`。
- 输出 `data/fig3_data.json`。

`mathematica/generate_fig3_single_case.wls`

- 历史单例：`a=0, Momega=0.6`。
- 角度网格是 2 度步长，不适合当前 281 点 paper-grid 图。

`mathematica/generate_fig3_w2_a0_fast.wls`

- 历史快速测试脚本。
- 不是当前主路径。

## 三、Python 求解/重求和文件

### `src/physics_fig3.py`

Fig.3 Python 数值核心。

主要职责：

- 读取 `fig3_data.json`、patch JSON、`fig3_sr_spherical.json`。
- 清理 Mathematica 高精度字符串。
- 计算 plain modal sum。
- 计算 Cesaro `alpha=2/5`。
- 计算 series reduction `k=2`。

当前实现细节：

- 若 `fig3_sr_spherical.json` 中有匹配 `(a, omega, k, theta_grid)` 的
  `f_sr(theta)`，优先使用这个 Appendix-C aligned `f`。
- `g(theta)` 当前仍用 modal sum。
- 若没有 spherical SR row，会 fallback 到 legacy spheroidal-index SR，并打印：
  `Falling back to legacy Fig3 SR ...`

### `main_fig3_paper_methods_grid.py`

当前推荐 Fig.3 画图入口。

输出：

- `results/Fig3_paper_methods_grid.png`
- `results/fig3_paper_methods_grid_metrics.json`

图像设置：

- 2x2 面板。
- 左列 `Momega=0.6`，右列 `Momega=2.0`。
- 上排 `a=0`，下排 `a=0.99M`。
- 红色实线：series reduction `k=2`。
- 蓝色点划线：Cesaro `alpha=2`。
- 黑色点划线：Cesaro `alpha=5`。
- `xlim = 30..180`，主刻度 `40..180`。
- `ylim = 0..80`，主刻度 `0,20,40,60,80`。

命令：

```powershell
python reproduce.py fig3-methods
```

### `scripts/audit_fig3_numeric_status.py`

当前推荐 Fig.3 数值审计入口。

输出：

- `results/fig3_numeric_audit.json`

审计内容：

- raw mode coverage。
- effective 281 点覆盖。
- spherical SR cache 覆盖。
- 每个 panel 的 raw min/max/p95。
- 有多少角度点超过 `y=80` 被图像裁剪。
- high-frequency precomputed SR 中 `f` 与 `g` 的拆分诊断。

命令：

```powershell
python reproduce.py fig3-audit
```

### 其他 Python 文件

`main_fig3.py`

- 旧 Fig.3 入口。
- 默认还是 `Momega=0.5,2.0` 的 1x2 历史复现路径。
- 不建议作为当前 2x2 paper-grid 主图入口。

`scripts/diag_fig3_w2_methods.py`

- 历史 `Momega=2` 方法比较诊断脚本。
- 当前主要保留作参考。

`reproduce.py`

当前 Fig.3 相关 target：

- `fig3`：旧 `main_fig3.py`
- `fig3-methods`：新 2x2 三方法图
- `fig3-audit`：数值覆盖审计

## 四、当前数据文件状态

### `data/fig3_data.json`

旧主缓存。

与当前 Fig.3 相关的覆盖：

- `a=0, Momega=0.5`: valid `l=2..74`
- `a=0, Momega=2.0`: valid 71 行，缺 `l=21,22`
- `a=0.99, Momega=0.5`: valid `l=2..51`
- `a=0.99, Momega=2.0`: valid `l=2..29`

### `data/fig3_patch_w2_a0_a099.json`

当前 patch 主缓存。

当前覆盖：

- `a=0, Momega=0.6`: 89 行，valid/281 点 `l=2..90`
- `a=0, Momega=2.0`: 13 行，valid/281 点 `l=2..12,21,22`
- `a=0.99, Momega=0.6`: 40 行，valid/281 点 `l=2..41`
- `a=0.99, Momega=2.0`: 73 行 valid；其中 281 点 `l=2..29`，71 点短网格 `l=30..74`

### `data/fig3_sr_spherical.json`

当前 spherical SR 缓存：

- `a=0, Momega=2.0, k=2`
- `a=0.99, Momega=2.0, k=2`

两个 row 都是旧缓存，`l_mode_max=60`。

## 五、当前输出

### `results/Fig3_paper_methods_grid.png`

当前最新 2x2 三方法图。

### `results/fig3_paper_methods_grid_metrics.json`

当前图的 raw 数值摘要：

- `Momega=0.6,a=0`: SR max `121.14`，`points_gt_80=9`
- `Momega=2.0,a=0`: SR max `265.23`，`points_gt_80=25`
- `Momega=0.6,a=0.99`: SR max `151.41`，`points_gt_80=27`
- `Momega=2.0,a=0.99`: SR max `121.55`，`points_gt_80=14`

### `results/fig3_numeric_audit.json`

当前完整审计结果。

## 六、当前主要瓶颈与风险

1. `Momega=0.6` panel 仍 fallback 到 legacy SR。
   原因是 `data/fig3_sr_spherical.json` 还没有 `Momega=0.6` 的 spherical-basis SR row。

2. `a=0.99, Momega=0.6` 高 `l` radial 计算很慢。
   当前已稳定到 `l=2..41`，继续 `l>=42` 时本地 Wolfram/BHPToolkit 长时间不落盘。

3. `a=0.99, Momega=2` 高 `l` angular densify 很慢或失败。
   当前 281 点覆盖到 `l=2..29`；`l=30..74` 保留旧 71 点 radial 数据。

4. `generate_fig3_sr_spherical.wls` 全 combo 重建本地超过 30 分钟未落盘。
   这是把 Fig.3 从“当前可复现图”推进到“论文级数值闭合图”的关键瓶颈。

5. 当前右列高频图的红色 SR 曲线仍有明显尖峰。
   这主要是数值覆盖/SR 缓存未完全闭合，不是画图坐标轴设置问题。

## 七、推荐后续顺序

1. 不再大批量跑 `a=0,Momega=2` 的已覆盖模态；
   旧数据已经覆盖大部分，只需保留补上的 `l=21,22`。

2. 专攻 `a=0.99,Momega=0.6,l>=42`：

```powershell
$env:FIG3_PATCH_A_LIST = "0.99"
$env:FIG3_PATCH_OMEGA_LIST = "0.6"
$env:FIG3_PATCH_LMIN = "42"
$env:FIG3_PATCH_LMAX = "42"
$env:FIG3_PATCH_BATCH = "1"
wolframscript -file mathematica/generate_fig3_patch_w2_highspin.wls
```

3. 分段 densify `a=0.99,Momega=2,l>=30`，必要时进一步调
   `NumTerms` 或改用更稳的 angular projection 策略。

4. 等 mode coverage 稳定后，再重建：

```powershell
$env:FIG3_SR_COMBOS = "0,0.6;0,2;0.99,0.6;0.99,2"
$env:FIG3_SR_LMAX = "88"
wolframscript -file mathematica/generate_fig3_sr_spherical.wls
python reproduce.py fig3-audit
python reproduce.py fig3-methods
```

