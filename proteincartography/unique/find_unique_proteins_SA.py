import pandas as pd

# ============================================================
# INPUT
# ============================================================

INPUT_FILE = "/Users/amcshan3/Desktop/Manuscripts/PLIP_Dpocket_Lipid_Puri_2025/Puri-et-al_Dpocket_PLIP_ProteinCartography_Lipid-Protein_Analysis/proteincartography/final_results_cluster-mode_saccharolipids/cluster-mode-saccharolipids_aggregated_features_pca_tsne.tsv"

df = pd.read_csv(INPUT_FILE, sep="\t")

# ============================================================
# CLEAN
# ============================================================

df["LeidenCluster"] = df["LeidenCluster"].fillna("NA")

# ============================================================
# FULL CLUSTER DISTRIBUTION
# ============================================================

lc_counts_full = df["LeidenCluster"].value_counts()
total = len(df)

# ============================================================
# METRICS
# ============================================================

top10_fraction = lc_counts_full.head(10).sum() / total * 100
top20_fraction = lc_counts_full.head(20).sum() / total * 100
top40_fraction = lc_counts_full.head(40).sum() / total * 100

print("\n==================== REDUNDANCY SUMMARY ====================\n")
print(f"Total structures: {total}\n")

print("LeidenCluster:")
print(f"  Unique clusters: {len(lc_counts_full)}")
print(f"  Median cluster size: {lc_counts_full.median():.2f}")
print(f"  Max cluster size: {lc_counts_full.max()}")
print(f"  Top 10 clusters cover: {top10_fraction:.2f}%")
print(f"  Top 20 clusters cover: {top20_fraction:.2f}%")
print(f"  Top 40 clusters cover: {top40_fraction:.2f}%\n")

# ============================================================
# DATA FOR PLOTTING
# ============================================================

top_n = 40
lc_counts = lc_counts_full.head(top_n)

lc_cum = lc_counts_full.cumsum() / total
lc_cum = lc_cum.head(top_n)

# ============================================================
# SVG EXPORT
# ============================================================

out_file = "LeidenCluster_representation_SA.svg"

width = 1200
height = 600

left_margin = 70
right_margin = 1120
bottom_margin = 550
top_margin = 50

plot_height = 450

bar_width = (right_margin - left_margin) / len(lc_counts)
max_val = lc_counts.max()

# ============================================================
# RIGHT AXIS (CUMULATIVE FRACTION 0–1)
# ============================================================

n_ticks = 6
right_tick_vals = [i / (n_ticks - 1) for i in range(n_ticks)]
right_tick_pos = [
    bottom_margin - (v * plot_height) for v in right_tick_vals
]

with open(out_file, "w") as f:
    f.write(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{height}">\n'
    )

    # --------------------------------------------------------
    # TITLE
    # --------------------------------------------------------

    f.write(
        '<text x="330" y="30" font-size="18" font-family="Arial">'
        'LeidenCluster Representation (Top 40)</text>\n'
    )

    # --------------------------------------------------------
    # AXES (LEFT + RIGHT)
    # --------------------------------------------------------

    # bottom axis
    f.write(
        f'<line x1="{left_margin}" y1="{bottom_margin}" '
        f'x2="{right_margin}" y2="{bottom_margin}" stroke="black"/>\n'
    )

    # left axis (counts)
    f.write(
        f'<line x1="{left_margin}" y1="{bottom_margin}" '
        f'x2="{left_margin}" y2="{top_margin}" stroke="black"/>\n'
    )

    # right axis (fraction)
    f.write(
        f'<line x1="{right_margin}" y1="{bottom_margin}" '
        f'x2="{right_margin}" y2="{top_margin}" stroke="black"/>\n'
    )

    # --------------------------------------------------------
    # LEFT Y-AXIS (COUNTS)
    # --------------------------------------------------------

    n_ticks = 6
    left_tick_vals = [int(i * max_val / (n_ticks - 1)) for i in range(n_ticks)]
    left_tick_pos = [
        bottom_margin - (v / max_val) * plot_height for v in left_tick_vals
    ]

    for tv, ty in zip(left_tick_vals, left_tick_pos):
        f.write(
            f'<line x1="{left_margin-5}" y1="{ty}" '
            f'x2="{left_margin}" y2="{ty}" stroke="black"/>\n'
        )
        f.write(
            f'<text x="{left_margin-10}" y="{ty+4}" font-size="10" '
            f'font-family="Arial" text-anchor="end">{tv}</text>\n'
        )

    # label
    f.write(
        f'<text x="15" y="{(bottom_margin+top_margin)/2}" font-size="14" '
        f'font-family="Arial" text-anchor="middle" dominant-baseline="middle" '
        f'transform="rotate(-90 15,{(bottom_margin+top_margin)/2})">'
        'Number of structures</text>\n'
    )

    # --------------------------------------------------------
    # RIGHT Y-AXIS (CUMULATIVE FRACTION)
    # --------------------------------------------------------

    for rv, ry in zip(right_tick_vals, right_tick_pos):
        f.write(
            f'<line x1="{right_margin}" y1="{ry}" '
            f'x2="{right_margin+5}" y2="{ry}" stroke="black"/>\n'
        )

        f.write(
            f'<text x="{right_margin+10}" y="{ry+4}" font-size="10" '
            f'font-family="Arial" text-anchor="start">{rv:.1f}</text>\n'
        )

    # label
    f.write(
        f'<text x="{width-15}" y="{(bottom_margin+top_margin)/2}" font-size="14" '
        f'font-family="Arial" text-anchor="middle" dominant-baseline="middle" '
        f'transform="rotate(90 {width-15},{(bottom_margin+top_margin)/2})">'
        'Cumulative fraction</text>\n'
    )

    # --------------------------------------------------------
    # BARS + X LABELS
    # --------------------------------------------------------

    for i, (cluster, val) in enumerate(lc_counts.items()):
        x = left_margin + i * bar_width
        bar_h = (val / max_val) * plot_height
        y = bottom_margin - bar_h

        f.write(
            f'<rect x="{x}" y="{y}" width="{bar_width*0.8}" height="{bar_h}" '
            f'style="fill:steelblue;stroke:black;stroke-width:1"/>\n'
        )

        f.write(
            f'<text x="{x}" y="565" font-size="10" font-family="Arial" '
            f'text-anchor="end" dominant-baseline="middle" '
            f'transform="rotate(-45 {x},565)">{cluster}</text>\n'
        )

    # --------------------------------------------------------
    # CUMULATIVE LINE (RIGHT AXIS SCALE)
    # --------------------------------------------------------

    points = []

    for i in range(len(lc_counts)):
        x = left_margin + i * bar_width + (bar_width * 0.4)
        y = bottom_margin - (lc_cum.iloc[i] * plot_height)
        points.append(f"{x},{y}")

    f.write(
        f'<polyline points="{" ".join(points)}" '
        f'style="fill:none;stroke:red;stroke-width:2"/>\n'
    )

    for i in range(len(lc_counts)):
        x = left_margin + i * bar_width + (bar_width * 0.4)
        y = bottom_margin - (lc_cum.iloc[i] * plot_height)

        f.write(
            f'<circle cx="{x}" cy="{y}" r="3" fill="red"/>\n'
        )

    # --------------------------------------------------------
    # CLOSE
    # --------------------------------------------------------

    f.write("</svg>")

print(f"\nSaved SVG: {out_file}")
print(f"Top 10 clusters cover {top10_fraction:.2f}%")
print(f"Top 20 clusters cover {top20_fraction:.2f}%")
print(f"Top 40 clusters cover {top40_fraction:.2f}%")
