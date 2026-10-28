#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATM alterations in the MSK-MET lung cancer cohort (Figure 7 and Results section
"ATM alterations in lung adenocarcinoma: frequency, genomic context and outcome
in the MSK-MET cohort"; Riedl, Reiss et al., Cancer Cell International).

Input
    data/ATM_MSK-MET_lung_extract.csv
        One row per LUAD/LUSC sample of the MSK-MET cohort (Nguyen et al., Cell 2022;
        cBioPortal study msk_met_2021, accessed 26 Oct 2026), one sample per patient.
        See README.md for the column definitions and the variant classification rules.

Usage
    python msk_met_atm_analysis.py [--data data/ATM_MSK-MET_lung_extract.csv] [--out results]

Output (written to --out)
    cohort_statistics.txt      all numbers reported in the manuscript
    coalterations.csv          odds ratios, 95% CI and Fisher p-values (Figure 7C)
    Figure_7_170mm.pdf / .png  Figure 7 (170 mm wide, 300 dpi for the PNG)
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

C_ANY, C_TR, C_WT = "#D95F02", "#7570B3", "#8C8C8C"

FEATURES = [("KRAS mutation", "KRAS_mut"), ("STK11 mutation", "STK11_mut"), ("KEAP1 mutation", "KEAP1_mut"),
            ("TP53 mutation", "TP53_mut"), ("EGFR mutation", "EGFR_mut"), ("BRAF mutation", "BRAF_mut"),
            ("MET exon 14", "MET_mut"), ("PIK3CA mutation", "PIK3CA_mut"), ("MET amplification", "MET_amp"),
            ("CDKN2A deletion", "CDKN2A_del")]


def odds_ratio_ci(a, b, c, d):
    """Odds ratio (a*d)/(b*c) with Woolf 95% CI; Haldane correction if any cell is 0."""
    t = np.array([[a, b], [c, d]], dtype=float)
    if (t == 0).any():
        t = t + 0.5
    orr = (t[0, 0] * t[1, 1]) / (t[0, 1] * t[1, 0])
    se = np.sqrt((1 / t).sum())
    return orr, np.exp(np.log(orr) - 1.96 * se), np.exp(np.log(orr) + 1.96 * se)


def analyse(df, out):
    lines = []
    def log(s=""):
        lines.append(s); print(s)

    luad = df[df.ONCOTREE_CODE == "LUAD"].copy()
    lusc = df[df.ONCOTREE_CODE == "LUSC"].copy()
    log("ATM alterations by histology")
    for name, g in [("LUAD", luad), ("LUSC", lusc)]:
        log(f"  {name}: n = {len(g):,}; ATM-altered (any) {int(g.ATM_any.sum())} ({100*g.ATM_any.mean():.1f}%); "
            f"ATM-truncated {int(g.ATM_truncating.sum())} ({100*g.ATM_truncating.mean():.1f}%)")
    for col, lab in [("ATM_any", "any alteration"), ("ATM_truncating", "truncating")]:
        p = fisher_exact([[luad[col].sum(), len(luad) - luad[col].sum()], [lusc[col].sum(), len(lusc) - lusc[col].sum()]])[1]
        log(f"  LUAD vs LUSC, {lab}: Fisher p = {p:.4f}")

    log("\nLUAD by sample type / site (any / truncating)")
    for s, g in luad.groupby("SAMPLE_TYPE"):
        log(f"  {s}: n = {len(g):,}; any {100*g.ATM_any.mean():.1f}%; truncating {100*g.ATM_truncating.mean():.1f}%")
    prim, met = luad[luad.SAMPLE_TYPE == "Primary"], luad[luad.SAMPLE_TYPE == "Metastasis"]
    for col, lab in [("ATM_any", "any"), ("ATM_truncating", "truncating")]:
        p = fisher_exact([[prim[col].sum(), len(prim) - prim[col].sum()], [met[col].sum(), len(met) - met[col].sum()]])[1]
        log(f"  primary vs metastasis, {lab}: Fisher p = {p:.2f}")
    for s, g in luad.groupby("SITE_GROUP"):
        log(f"  {s}: n = {len(g)}; any {100*g.ATM_any.mean():.1f}%; truncating {100*g.ATM_truncating.mean():.1f}%")

    log("\nCo-alterations in LUAD (vs. ATM wild-type)")
    wt = luad[luad.ATM_any == 0]; anyg = luad[luad.ATM_any == 1]; trg = luad[luad.ATM_truncating == 1]
    log(f"  ATM wild-type n = {len(wt):,}; ATM-altered n = {len(anyg)}; ATM-truncated n = {len(trg)}")
    rows = []
    for name, col in FEATURES:
        b = wt[col].fillna(0).astype(int)
        for lab, g in [("ATM-truncated", trg), ("ATM-altered", anyg)]:
            a = g[col].fillna(0).astype(int)
            orr, lo, hi = odds_ratio_ci(a.sum(), len(a) - a.sum(), b.sum(), len(b) - b.sum())
            p = fisher_exact([[a.sum(), len(a) - a.sum()], [b.sum(), len(b) - b.sum()]])[1]
            rows.append(dict(feature=name, group=lab, pct_group=100 * a.mean(), pct_wt=100 * b.mean(),
                             odds_ratio=orr, ci_low=lo, ci_high=hi, fisher_p=p))
            log(f"  {name:18s} {lab:14s} {100*a.mean():5.1f}% vs {100*b.mean():5.1f}%  OR {orr:.2f} ({lo:.2f}-{hi:.2f})  p = {p:.1e}")
    co = pd.DataFrame(rows)
    co.to_csv(os.path.join(out, "coalterations.csv"), index=False)

    log("\nTumour mutational burden (non-synonymous mutations per Mb)")
    p = mannwhitneyu(anyg.TMB_NONSYNONYMOUS.dropna(), wt.TMB_NONSYNONYMOUS.dropna())[1]
    log(f"  median: wild-type {wt.TMB_NONSYNONYMOUS.median():.1f}; ATM-altered {anyg.TMB_NONSYNONYMOUS.median():.1f}; "
        f"ATM-truncated {trg.TMB_NONSYNONYMOUS.median():.1f}; Mann-Whitney U (altered vs wild-type) p = {p:.1e}")

    log("\nOverall survival from sequencing (LUAD)")
    s = luad.dropna(subset=["OS_MONTHS", "OS_STATUS"]).copy()
    s["event"] = (s.OS_STATUS == "1:DECEASED").astype(int)
    s["grp"] = np.where(s.ATM_truncating == 1, "truncated", np.where(s.ATM_any == 1, "missense", "wt"))
    km = KaplanMeierFitter(); med = {}
    for g in ["wt", "missense", "truncated"]:
        gg = s[s.grp == g]; km.fit(gg.OS_MONTHS, gg.event); med[g] = km.median_survival_time_
        log(f"  {g}: n = {len(gg):,}; events = {int(gg.event.sum())}; median OS = {med[g]:.1f} months")
    wtg = s[s.grp == "wt"]
    p_any = logrank_test(s[s.ATM_any == 1].OS_MONTHS, wtg.OS_MONTHS, s[s.ATM_any == 1].event, wtg.event).p_value
    p_tr = logrank_test(s[s.grp == "truncated"].OS_MONTHS, wtg.OS_MONTHS, s[s.grp == "truncated"].event, wtg.event).p_value
    log(f"  log-rank vs wild-type: any alteration p = {p_any:.3f}; truncated p = {p_tr:.3f}; n with follow-up = {len(s):,}")
    with open(os.path.join(out, "cohort_statistics.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return luad, lusc, co, s, med, p_any, p_tr


def figure7(luad, lusc, co, s, med, p_any, p_tr, out):
    import matplotlib.font_manager as fm
    avail = {f.name for f in fm.fontManager.ttflist}
    family = next((f for f in ["Arial", "Liberation Sans", "DejaVu Sans"] if f in avail), "sans-serif")
    plt.rcParams.update({"font.family": family, "font.size": 7,
                         "axes.linewidth": 1.0, "axes.labelweight": "bold", "axes.titleweight": "bold",
                         "xtick.major.width": 1.0, "ytick.major.width": 1.0, "xtick.direction": "out",
                         "ytick.direction": "out", "xtick.major.size": 3.5, "ytick.major.size": 3.5,
                         "legend.fontsize": 6.5, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(6.69, 6.6))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 1.15], width_ratios=[1.0, 1.35], hspace=0.5, wspace=0.42,
                           left=0.09, right=0.98, top=0.95, bottom=0.08)
    # A: frequency
    ax = fig.add_subplot(gs[0, 0]); w = 0.36
    for i, (lab, g) in enumerate([("Lung adeno-\ncarcinoma\n(n = 4,064)", luad), ("Squamous cell\ncarcinoma\n(n = 537)", lusc)]):
        pa, pt = 100 * g.ATM_any.mean(), 100 * g.ATM_truncating.mean()
        ax.bar(i - w / 2, pa, w, color=C_ANY, edgecolor="black", linewidth=0.8, label="ATM-altered (any)" if i == 0 else None)
        ax.bar(i + w / 2, pt, w, color=C_TR, edgecolor="black", linewidth=0.8, label="ATM-truncated" if i == 0 else None)
        ax.text(i - w / 2, pa + 0.25, f"{pa:.1f}%\n(n = {int(g.ATM_any.sum())})", ha="center", va="bottom", fontsize=6)
        ax.text(i + w / 2, pt + 0.25, f"{pt:.1f}%\n(n = {int(g.ATM_truncating.sum())})", ha="center", va="bottom", fontsize=6)
    pa = fisher_exact([[luad.ATM_any.sum(), len(luad) - luad.ATM_any.sum()], [lusc.ATM_any.sum(), len(lusc) - lusc.ATM_any.sum()]])[1]
    pt = fisher_exact([[luad.ATM_truncating.sum(), len(luad) - luad.ATM_truncating.sum()], [lusc.ATM_truncating.sum(), len(lusc) - lusc.ATM_truncating.sum()]])[1]
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Lung adeno-\ncarcinoma\n(n = 4,064)", "Squamous cell\ncarcinoma\n(n = 537)"], fontsize=6.5)
    ax.set_ylabel("Samples with ATM alteration (%)"); ax.set_ylim(0, 11); ax.set_title("ATM alterations by histology", fontsize=7.5)
    ax.text(0.02, 0.98, f"LUAD vs. LUSC (Fisher's exact test):\nany p = {pa:.3f}; truncated p = {pt:.4f}", transform=ax.transAxes, va="top", fontsize=6)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.80), frameon=False, fontsize=6, handlelength=1.0)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    ax.text(-0.28, 1.05, "A", transform=ax.transAxes, fontsize=13, fontweight="bold")
    # B: survival
    ax = fig.add_subplot(gs[0, 1]); km = KaplanMeierFitter()
    for g, col, lab in [("wt", C_WT, "ATM wild-type"), ("missense", C_ANY, "ATM-altered, missense/in-frame only"), ("truncated", C_TR, "ATM-truncated")]:
        gg = s[s.grp == g]; km.fit(gg.OS_MONTHS, gg.event, label=f"{lab} (n = {len(gg):,})")
        km.plot_survival_function(ax=ax, ci_show=False, color=col, lw=1.6, show_censors=False)
    ax.set_xlim(0, 72); ax.set_ylim(0, 1.0); ax.set_xlabel("Time from sequencing (months)"); ax.set_ylabel("Overall survival probability")
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    ax.text(0.02, 0.06, f"log-rank vs. wild-type: any alteration p = {p_any:.2f}; truncated p = {p_tr:.2f}\n"
            f"median OS: wild-type {med['wt']:.1f}, missense/in-frame {med['missense']:.1f}, truncated {med['truncated']:.1f} months",
            transform=ax.transAxes, fontsize=5.8, va="bottom")
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    ax.set_title(f"Overall survival by ATM status (LUAD, n = {len(s):,})", fontsize=7.5)
    ax.text(-0.2, 1.05, "B", transform=ax.transAxes, fontsize=13, fontweight="bold")
    # C: co-alterations
    ax = fig.add_subplot(gs[1, :])
    names = [f for f, _ in FEATURES]; y = np.arange(len(names))[::-1]
    for i, name in enumerate(names):
        for lab, col, off, mk in [("ATM-truncated", C_TR, 0.18, "o"), ("ATM-altered", C_ANY, -0.18, "s")]:
            r = co[(co.feature == name) & (co.group == lab)].iloc[0]
            ax.plot([r.ci_low, r.ci_high], [y[i] + off, y[i] + off], color=col, lw=1.4)
            ax.plot(r.odds_ratio, y[i] + off, mk, color=col, ms=5, mec="black", mew=0.5)
    ax.axvline(1, color="k", ls="--", lw=0.7); ax.set_xscale("log"); ax.set_xlim(0.03, 12)
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4, 8]); ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4", "8"])
    ax.set_yticks(y); ax.set_yticklabels(names); ax.set_ylim(-0.7, len(names) - 0.3)
    ax.set_xlabel("Odds ratio vs. ATM wild-type LUAD (95% CI)")
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    n_wt = int((luad.ATM_any == 0).sum()); n_any = int(luad.ATM_any.sum()); n_tr = int(luad.ATM_truncating.sum())
    ax.plot([], [], "o", color=C_TR, ms=4, label=f"ATM-truncated (nonsense, frameshift, splice; n = {n_tr}) vs. ATM wild-type (n = {n_wt:,})")
    ax.plot([], [], "s", color=C_ANY, ms=4, label=f"ATM-altered (any non-synonymous variant; n = {n_any}) vs. ATM wild-type")
    ax.legend(frameon=False, fontsize=6, loc="lower left", bbox_to_anchor=(0.0, 1.0))
    xt = 13.5
    for x, t in [(xt, "group vs. WT"), (xt * 2.6, "OR"), (xt * 3.6, "Fisher p")]:
        ax.text(x, len(names) - 0.25, t, fontsize=6, fontweight="bold", va="center", clip_on=False)
    fp = lambda p: f"{p:.2f}" if p >= 0.01 else f"{p:.0e}".replace("e-0", "e-")
    for i, name in enumerate(names):
        for lab, col, off in [("ATM-truncated", C_TR, 0.18), ("ATM-altered", C_ANY, -0.18)]:
            r = co[(co.feature == name) & (co.group == lab)].iloc[0]
            ax.text(xt, y[i] + off, f"{r.pct_group:.1f}% vs {r.pct_wt:.1f}%", fontsize=5.6, color=col, va="center", clip_on=False)
            ax.text(xt * 2.6, y[i] + off, f"{r.odds_ratio:.2f}", fontsize=5.6, color=col, va="center", clip_on=False)
            ax.text(xt * 3.6, y[i] + off, fp(r.fisher_p), fontsize=5.6, color=col, va="center", clip_on=False)
    ax.text(-0.22, 1.05, "C", transform=ax.transAxes, fontsize=13, fontweight="bold")
    pos = ax.get_position(); ax.set_position([pos.x0 + 0.06, pos.y0, pos.width * 0.62, pos.height])
    fig.savefig(os.path.join(out, "Figure_7_170mm.pdf")); fig.savefig(os.path.join(out, "Figure_7_170mm.png"), dpi=300)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="data/ATM_MSK-MET_lung_extract.csv")
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    df = pd.read_csv(a.data)
    luad, lusc, co, s, med, p_any, p_tr = analyse(df, a.out)
    figure7(luad, lusc, co, s, med, p_any, p_tr, a.out)
    print(f"\nFigure 7 and statistics written to {a.out}/")


if __name__ == "__main__":
    main()
