"""Generate assets/iron-yield-ratios.svg directly from frozen evidence JSONs."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png", type=Path, default=None, help="Optional PNG preview output path")
    args = parser.parse_args()

    repo = Path(__file__).resolve().parent.parent
    out_svg = repo / "assets" / "iron-yield-ratios.svg"

    plan = json.loads((repo / "evidence/records/formal-transport-plan.json").read_text(encoding="utf-8"))
    fe_m01 = json.loads((repo / "evidence/units/Fe-M01.json").read_text(encoding="utf-8"))
    fe_c01 = json.loads((repo / "evidence/units/Fe-C01.json").read_text(encoding="utf-8"))

    ratios_by_pair = {
        (r["numerator"], r["denominator"]): r for r in plan["reference_table"]["ratios"]
    }
    r_def = ratios_by_pair[("g4-default", "fluka-default")]
    r_on = ratios_by_pair[("g4-default", "fluka-coalesce")]
    r_off = ratios_by_pair[("g4-cluster-off", "fluka-default")]
    r_rev = ratios_by_pair[("g4-cluster-off", "fluka-coalesce")]

    ref_r0 = plan["selected_reference"]["ratio"]
    ref_2se = 2.0 * plan["selected_reference"]["mc_se"]

    rows = [
        (6.2, "Defaults: G4 on, FLUKA off", r_def["ratio"], 2 * r_def["mc_se"], "{:.3f} ± {:.3f}", "#64748B", "o", False),
        (5.2, "Both on: preselected ref.", r_on["ratio"], 2 * r_on["mc_se"], "{:.3f} ± {:.3f}", "#1D4ED8", "D", True),
        (4.2, "Both off: G4 off, FLUKA off", r_off["ratio"], 2 * r_off["mc_se"], "{:.3f} ± {:.3f}", "#64748B", "o", False),
        (3.2, "Reversed: G4 off, FLUKA on", r_rev["ratio"], 2 * r_rev["mc_se"], "{:.3f} ± {:.3f}", "#64748B", "o", False),
        (1.7, "Fe-M01: both on at init", fe_m01["ratio"], 2 * fe_m01["mc_se"], "{:.4f} ± {:.4f}", "#C2410C", "s", True),
        (0.7, "Fe-C01: both on at init", fe_c01["ratio"], 2 * fe_c01["mc_se"], "{:.4f} ± {:.4f}", "#C2410C", "s", True),
    ]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 10,
        "svg.fonttype": "none",
        "svg.hashsalt": "mcpair-iron-ratios-v1",
        "axes.edgecolor": "#475569",
        "axes.linewidth": 0.9,
    })

    fig, ax = plt.subplots(figsize=(8.8, 4.4), dpi=180)
    fig.patch.set_facecolor("#FFFFFF")
    ax.set_facecolor("#FFFFFF")

    ax.axvspan(ref_r0 - ref_2se, ref_r0 + ref_2se, color="#DBEAFE", alpha=0.65, zorder=1)
    ax.axvline(ref_r0, color="#1D4ED8", linestyle="--", linewidth=1.1, alpha=0.85, zorder=2)
    ax.axvline(1.0, color="#94A3B8", linestyle=":", linewidth=1.1, zorder=2)
    ax.axhline(2.55, color="#CBD5E1", linestyle="-", linewidth=0.9, zorder=1)

    hdr_box = dict(facecolor="#FFFFFF", edgecolor="none", pad=1.8, alpha=0.92)
    ax.text(0.935, 6.80, "Earlier reference runs", fontsize=9.2, fontweight="bold", color="#334155", va="bottom", ha="left", bbox=hdr_box, zorder=3)
    ax.text(0.935, 2.16, "Formal LLM-written pairs", fontsize=9.2, fontweight="bold", color="#9A3412", va="bottom", ha="left", bbox=hdr_box, zorder=3)

    for y, label, val, err, fmt, color, marker, is_bold in rows:
        ax.errorbar(
            val,
            y,
            xerr=err,
            fmt=marker,
            color=color,
            ecolor=color,
            elinewidth=1.6,
            capsize=3.5,
            capthick=1.4,
            markersize=6.2 if marker != "D" else 6.0,
            zorder=4,
        )
        ax.text(
            val + err + 0.018,
            y,
            fmt.format(val, err),
            va="center",
            ha="left",
            fontsize=9.3,
            fontweight="bold" if is_bold else "normal",
            color="#0F172A" if is_bold else "#334155",
            zorder=5,
        )

    ax.text(1.0, 7.26, "1.0", ha="center", va="bottom", fontsize=8.8, color="#64748B")
    ax.text(ref_r0, 7.26, "Reference ±2 MCSE", ha="center", va="bottom", fontsize=8.8, fontweight="bold", color="#1D4ED8")

    ax.set_yticks([r[0] for r in rows])
    ax.set_yticklabels([r[1] for r in rows])
    for tick_label, row in zip(ax.get_yticklabels(), rows):
        tick_label.set_fontsize(9.6)
        if row[7]:
            tick_label.set_fontweight("bold")
            tick_label.set_color("#0F172A")
        else:
            tick_label.set_color("#334155")

    ax.set_xlim(0.92, 2.04)
    ax.set_ylim(0.15, 7.22)
    ax.set_xlabel("Geant4 / FLUKA free-neutron yield ratio (±2 MCSE)", fontsize=10, color="#0F172A", labelpad=7)

    ax.grid(axis="x", color="#E2E8F0", linestyle="-", linewidth=0.7, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.tick_params(axis="x", colors="#334155")

    fig.suptitle(
        "Thin natural iron: 800 MeV protons, 4π, 35.18–107.29 MeV neutron kinetic energy",
        x=0.025,
        y=0.965,
        ha="left",
        fontsize=10.8,
        fontweight="bold",
        color="#0F172A",
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90))
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_svg, format="svg", facecolor=fig.get_facecolor(), metadata={"Date": None, "Creator": "MCPair plot_iron_ratios.py"})
    print(f"Wrote {out_svg} ({out_svg.stat().st_size} bytes)")
    if args.png is not None:
        args.png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.png, format="png", dpi=180, facecolor=fig.get_facecolor())
        print(f"Wrote {args.png} ({args.png.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
