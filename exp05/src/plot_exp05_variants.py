"""Two traceable renderings of the theory-derived root-locus sketch.

Uses existing Octave selection/trajectory CSVs and an independent algebraic
branch equation. The pen rendering is a program-drawn style, not a scan.
Run with the existing Anaconda Python; no added dependencies are required.
"""
from pathlib import Path
import json
import os

ROOT = Path(__file__).resolve().parents[2]
PRIVATE = ROOT / ".local" / "exp05"
os.environ["MPLCONFIGDIR"] = str(PRIVATE / "mplconfig")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = PRIVATE / "results"
plt.rcParams.update({
    "font.family": ["Microsoft YaHei", "DejaVu Sans"],
    "font.size": 11, "axes.unicode_minus": False,
    "svg.fonttype": "none", "axes.spines.top": False,
    "axes.spines.right": False,
})


def draw(pen=False):
    selection = np.genfromtxt(OUT / "selection.csv", delimiter=",", names=True)
    software = np.loadtxt(OUT / "root_locus.csv", delimiter=",", skiprows=1)
    xb = -2 + 2 / np.sqrt(3)
    kb = -xb * (xb + 2) * (xb + 4)
    # Quadratic spacing resolves the initially vertical tangent at breakaway.
    t = np.linspace(0, 1, 1300)
    x = xb + (1.15 - xb) * t**2
    y = np.sqrt(np.maximum(0, 3*x*x + 12*x + 8))
    curve = "#254967" if pen else "#1269a5"
    fig = plt.figure(figsize=(12, 7.5), layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=[3.1, 1.15])
    ax = fig.add_subplot(gs[0]); notes = fig.add_subplot(gs[1]); notes.axis("off")
    ax.set(xlim=(-7.7, 2.0), ylim=(-4.9, 4.9),
           xlabel=r"实部 $\mathrm{Re}(s)$ / $\mathrm{s}^{-1}$",
           ylabel=r"虚部 $\mathrm{Im}(s)$ / $\mathrm{s}^{-1}$")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(np.arange(-7, 3)); ax.set_yticks(np.arange(-4, 5))
    ax.axhline(0, color="#555b61", lw=0.9, zorder=0)
    ax.axvline(0, color="#555b61", lw=0.9, zorder=0)
    ax.grid(color="#d9dfe4", linewidth=0.6, alpha=0.7)

    def line(xx, yy, **kw):
        obj, = ax.plot(xx, yy, **kw)
        if pen:
            # Only the strokes change visually; axes and exact point markers
            # remain stable. This does not perturb the numerical source data.
            obj.set_sketch_params(scale=0.65, length=110, randomness=1.4)
        return obj

    for angle in (np.pi/3, np.pi, 5*np.pi/3):
        u = np.linspace(0, 12, 200)
        line(-2+u*np.cos(angle), u*np.sin(angle), color="#93999f",
             lw=1.1, linestyle=(0, (5, 4)), zorder=1)
    # Positive-gain complex branches, independently derived from Im(q)=0.
    for sign in (-1, 1):
        line(x, sign*y, color=curve, lw=2.1, zorder=3)
    line(np.linspace(-7.7, -4, 300), np.zeros(300), color=curve, lw=2.1)
    line(np.linspace(-2, 0, 300), np.zeros(300), color=curve, lw=2.1)
    # Independent software output provides a restrained visual cross-check.
    if not pen:
        for j in range(3):
            ax.scatter(software[::5, 1+2*j], software[::5, 2+2*j],
                       s=8, color=curve, alpha=0.22, zorder=2)
    arrows = [((-5.1, 0), (-5.9, 0)), ((-.12, 0), (-.52, 0)),
              ((-1.88, 0), (-1.3, 0))]
    for sign in (-1, 1):
        arrows.append(((-.45, sign*np.sqrt(3*(-.45)**2+12*(-.45)+8)),
                       (-.22, sign*np.sqrt(3*(-.22)**2+12*(-.22)+8))))
    for a, b in arrows:
        ax.annotate("", xy=b, xytext=a,
                    arrowprops={"arrowstyle": "-|>", "color": curve,
                                "lw": 1.5, "mutation_scale": 14})
    ax.scatter([-4, -2, 0], [0, 0, 0], marker="x", s=95,
               linewidths=2.1, color="#17212b", zorder=6)
    ax.scatter([xb], [0], s=52, color="#ad6c1b", zorder=7)
    ax.scatter([0, 0], [np.sqrt(8), -np.sqrt(8)], s=65,
               edgecolors="#b74136", facecolors="white", linewidths=1.8, zorder=7)
    ax.scatter([-2], [0], marker="s", s=44, facecolors="white",
               edgecolors="#626971", linewidths=1.3, zorder=5)
    sx, sy = float(selection["selected_real"]), float(selection["selected_imag"])
    ax.scatter([sx, sx], [sy, -sy], s=53, color="#19846b", zorder=7)
    ax.annotate(r"$K=48,\ s=+j\sqrt{8}$", (0, np.sqrt(8)),
                xytext=(-3.7, 3.55), arrowprops={"arrowstyle": "-", "color": "#b74136"},
                color="#a6352c", fontsize=10)
    ax.annotate(r"$s_b=-0.84530$", (xb, 0), xytext=(-4.1, -.9),
                arrowprops={"arrowstyle": "-", "color": "#ad6c1b"}, color="#8e5715")
    ax.annotate(r"$\zeta=0.5$ 设计点", (sx, sy), xytext=(-4.5, 2.0),
                arrowprops={"arrowstyle": "-", "color": "#19846b"}, color="#16735d")
    note_font = ["KaiTi", "DejaVu Sans"] if pen else ["Microsoft YaHei", "DejaVu Sans"]
    notes.text(0, .97, "对象与关键计算", transform=notes.transAxes,
               va="top", fontsize=14, weight="bold", fontfamily=note_font)
    text = (
        "$L(s)=K/[s(s+2)(s+4)]$\n\n"
        "三条分支，无有限零点\n"
        "起点：0，−2，−4\n\n"
        "实轴区段\n"
        r"$(-\infty,-4),\ (-2,0)$" "\n\n"
        "渐近线交点：−2\n"
        "角度：60°、180°、300°\n\n"
        f"分离点增益：{kb:.5f}\n"
        "虚轴交点增益：48\n"
        r"$s=\pm j\sqrt{8}$" "\n\n"
        "稳定范围：0 < K < 48\n\n"
        r"$\zeta=0.5$ 设计点" "\n"
        r"$K=224/27\approx8.29630$" "\n"
        r"$s=-2/3\pm j\,2\sqrt{3}/3$"
    )
    notes.text(0, .89, text, transform=notes.transAxes, va="top",
               fontsize=10.7, linespacing=1.48, fontfamily=note_font)
    variant = "笔绘风格草图" if pen else "标准 Matplotlib 图"
    fig.suptitle(f"实验五｜根轨迹关键结构 · {variant}", fontsize=17)
    ax.legend(handles=[Line2D([], [], color=curve, lw=2, label="正增益根轨迹"),
                       Line2D([], [], color="#93999f", lw=1, ls="--", label="渐近线"),
                       Line2D([], [], color="#17212b", marker="x", ls="none", label="开环极点")],
              loc="lower left", fontsize=9, framealpha=.95)
    fig.supxlabel("依据特征方程推导并程序绘制；箭头表示 K 增大方向。", fontsize=10)
    stem = "A_root_locus_penstyle" if pen else "A_root_locus_matplotlib"
    fig.savefig(OUT / f"{stem}.png", dpi=220)
    fig.savefig(OUT / f"{stem}.svg")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    draw(False); draw(True)
    provenance = {
        "figure_stems": ["A_root_locus_matplotlib", "A_root_locus_penstyle"],
        "code": "exp05/src/plot_exp05_variants.py",
        "inputs": [".local/exp05/results/root_locus.csv", ".local/exp05/results/selection.csv"],
        "branch_equation": "Im(q)=0: y^2=3*x^2+12*x+8, x>=-2+2/sqrt(3)",
        "pen_style": "Program-drawn strokes; stable axes and exact key markers; not a paper scan",
    }
    (OUT / "sketch_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    print("Saved both root-locus variants (PNG/SVG) and provenance.")


if __name__ == "__main__":
    main()
