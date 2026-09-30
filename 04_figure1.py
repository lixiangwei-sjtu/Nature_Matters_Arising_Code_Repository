"""Reproduce the manuscript Figure 1 from non-disclosive locked estimates."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results", "locked_main_results.csv")
out = pd.read_csv(RESULTS)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = True
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

HRS_COLOR = "#0072B2"
NHANES_COLOR = "#D55E00"
NULL_COLOR = "0.80"

fig = plt.figure(figsize=(14.0, 7.2), facecolor="white")
PANEL_LEFT = 0.075
PANEL_WIDTH = 0.85
HEADER_OFFSET = 0.018
PANEL_GAP = 0.105
PANEL_A_HEIGHT = 0.145
PANEL_B_HEIGHT = 0.155
PANEL_C_HEIGHT = 0.245
PANEL_C_BOTTOM = 0.065
PANEL_B_BOTTOM = PANEL_C_BOTTOM + PANEL_C_HEIGHT + PANEL_GAP
PANEL_A_BOTTOM = PANEL_B_BOTTOM + PANEL_B_HEIGHT + PANEL_GAP
COL_LABEL, COL_FOREST, COL_EST, COL_P = 0.245, 0.320, 0.335, 0.080
COL_GAP = 0.0025
XMIN, XMAX = -0.13, 0.080
XTICKS = [-0.12, -0.08, -0.04, 0.00, 0.04]

# panel a
ax_a = fig.add_axes([PANEL_LEFT, PANEL_A_BOTTOM, PANEL_WIDTH, PANEL_A_HEIGHT])
ax_a.set_xlim(0, 1); ax_a.set_ylim(0, 1); ax_a.axis("off")
fig.text(PANEL_LEFT, PANEL_A_BOTTOM + PANEL_A_HEIGHT + HEADER_OFFSET, "a", fontsize=11, fontweight="bold", va="center")
fig.text(PANEL_LEFT + 0.025, PANEL_A_BOTTOM + PANEL_A_HEIGHT + HEADER_OFFSET, "Drosophila mechanism", fontsize=10.5, fontweight="bold", va="center")
box_w, box_h, box_y = 0.17, 0.34, 0.55
centres = [0.09, 0.36, 0.63, 0.90]
labels = ["Larval Phe/Tyr\nrestriction", "Lsp2 ↓", "Early-adult\ntranslation ↓", "Lifespan ↑"]
for xc, lab in zip(centres, labels):
    ax_a.add_patch(patches.Rectangle((xc-box_w/2, box_y-box_h/2), box_w, box_h, linewidth=0.75, edgecolor="black", facecolor="white"))
    ax_a.text(xc, box_y, lab, ha="center", va="center", fontsize=8.8, linespacing=0.94)
for i in range(3):
    ax_a.annotate("", xy=(centres[i+1]-box_w/2, box_y), xytext=(centres[i]+box_w/2, box_y), arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black", mutation_scale=9))
ax_a.text(0.50, 0.10, "Developmental intervention", ha="center", va="center", fontsize=7.8, fontstyle="italic", color="0.25")


def draw_panel(letter, title, data, bottom, height, y_positions, ymin, ymax, header_y):
    left, width = PANEL_LEFT, PANEL_WIDTH
    fig.text(left, bottom+height+HEADER_OFFSET, letter, fontsize=11, fontweight="bold", va="center")
    fig.text(left+0.025, bottom+height+HEADER_OFFSET, title, fontsize=10.5, fontweight="bold", va="center")
    lw, fw, ew, pw = [width*x for x in (COL_LABEL, COL_FOREST, COL_EST, COL_P)]
    ax_lab = fig.add_axes([left, bottom, lw, height]); ax_lab.axis("off"); ax_lab.set_ylim(ymin,ymax); ax_lab.set_xlim(0,1)
    ax_for = fig.add_axes([left+lw+COL_GAP, bottom, fw, height]); ax_for.set_ylim(ymin,ymax); ax_for.set_xlim(XMIN,XMAX)
    ax_est = fig.add_axes([left+lw+COL_GAP+fw+COL_GAP, bottom, ew, height]); ax_est.axis("off"); ax_est.set_ylim(ymin,ymax); ax_est.set_xlim(0,1)
    ax_p = fig.add_axes([left+lw+COL_GAP+fw+COL_GAP+ew+COL_GAP, bottom, pw, height]); ax_p.axis("off"); ax_p.set_ylim(ymin,ymax); ax_p.set_xlim(0,1)
    ax_for.axvline(0, color=NULL_COLOR, linestyle=(0,(3,3)), linewidth=0.65)
    for j, (_, r) in enumerate(data.iterrows()):
        y = y_positions[j]
        label = f"{r.cohort}  {r.exposure}"
        ax_lab.text(0.98, y, label, ha="right", va="center", fontsize=8.4)
        color = HRS_COLOR if r.cohort == "HRS" else NHANES_COLOR
        mfc = color if r.cohort == "HRS" else "white"
        xerr = np.array([[r.beta-r.ci_low],[r.ci_high-r.beta]])
        ax_for.errorbar(r.beta, y, xerr=xerr, fmt="o", markersize=5.1, color=color, ecolor=color, markerfacecolor=mfc, markeredgecolor=color, markeredgewidth=1.2, elinewidth=1.15, capsize=3.4, capthick=1.0)
        sign = "+" if r.beta >= 0 else "−"
        b = abs(r.beta); lo = r.ci_low; hi = r.ci_high
        txt = f"{sign}{b:.4f} ({lo:+.4f}, {hi:+.4f})".replace("-","−")
        ax_est.text(0.04, y, txt, ha="left", va="center", fontsize=7.6)
        ax_p.text(0.48, y, f"{r.p_value:.3f}", ha="center", va="center", fontsize=7.8)
    ax_est.text(0.04, header_y, "β (95% CI)", ha="left", va="center", fontsize=8.0)
    ax_p.text(0.48, header_y, "P", ha="center", va="center", fontsize=8.0, fontstyle="italic")
    ax_for.set_xticks(XTICKS); ax_for.set_xticklabels([f"{x:.2f}" for x in XTICKS], fontsize=7.5); ax_for.set_yticks([])
    ax_for.set_xlabel(r"Standardized $\beta$", fontsize=8.8, labelpad=6)
    for s in ["top","right","left"]: ax_for.spines[s].set_visible(False)
    ax_for.spines["bottom"].set_linewidth(0.8)

b = out[out.panel=="b"].reset_index(drop=True)
c = out[out.panel=="c"].reset_index(drop=True)
draw_panel("b", "Direct Phe/Tyr test", b, PANEL_B_BOTTOM, PANEL_B_HEIGHT, [2,1], 0.45, 3.35, 3.00)
draw_panel("c", "Albumin-informed compositional test", c, PANEL_C_BOTTOM, PANEL_C_HEIGHT, [4,3,2,1], 0.45, 5.35, 5.00)

fig.savefig("Figure1_Nature_Matters_Arising.pdf", facecolor="white")
fig.savefig("Figure1_Nature_Matters_Arising.png", dpi=600, facecolor="white")
print("Saved Figure1_Nature_Matters_Arising.pdf and .png")
