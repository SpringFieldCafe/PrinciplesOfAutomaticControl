# 自动控制原理

本项目分为 [上机实验](computer-labs/) 和 [书面作业](written-assignments/) 两部分。上机实验整理 8 个实验的原始课件、数值分析与绘图源码；八个实验的代码均已实际运行，实验一另附三份可编辑的 Xcos 闭环模型。书面作业目录预留，待题目材料到位后补充。

| 目录 | 实验 | Octave 数值入口 |
| --- | --- | --- |
| [exp01](computer-labs/exp01/) | 控制系统仿真模型建立 | `run_exp01()` |
| [exp02](computer-labs/exp02/) | 反馈系统特性分析 | `run_exp02('compute')` |
| [exp03](computer-labs/exp03/) | 控制系统性能指标分析 | `run_exp03(false)` |
| [exp04](computer-labs/exp04/) | 控制系统稳定性分析 | `run_exp04()` |
| [exp05](computer-labs/exp05/) | 控制系统根轨迹法分析 | `run_exp05(false)` |
| [exp06](computer-labs/exp06/) | 控制系统频率特性分析 | `run_exp06('all')` |
| [exp07](computer-labs/exp07/) | 控制系统频域稳定性分析 | `run_exp07('compute')` |
| [exp08](computer-labs/exp08/) | 反馈控制系统设计 | `run_exp08('nominal')` |

数值分析使用 GNU Octave 11.3.0 与 control 4.2.1；实验一框图使用 Scilab 2026.1.0 / Xcos。部分图表使用 Python 3.12.7、NumPy 1.26.4、Matplotlib 3.9.2 绘制。

各实验的 `doc/` 保存原始课件，`src/` 保存代码及模型。运行前先阅读课件，并将目标实验的 `src/` 加入 Octave 搜索路径，例如在项目根目录执行 `addpath('computer-labs/exp01/src'); run_exp01();`。源码根据自身位置定位资源，结果输出至 `computer-labs/.local/expXX/results/`；实验六首次运行前需创建该结果目录。

实验二、七计算后调用同一入口的 `'plot'` 模式；实验三、五分别调用 `plot_exp03()`、`plot_exp05()`。独立 Python 绘图脚本位于相应 `src/` 中。实验八的 `'robust'` 模式计算参数扰动，实验六的 `verify_exp06()` 核验辨识结果；实验一的 `.zcos` 模型可在 Xcos 中打开、编辑和运行。

Word 报告、生成图表与数据、交付 ZIP、软件环境和本地协作材料保留在本地，由 Git 忽略规则排除，不随仓库发布。原始课件及其备注保持原样。

本仓库的代码、说明和实验 PPT 均采用 [MIT 许可证](LICENSE)。课件中引用的外部资料保留其原有权属。
