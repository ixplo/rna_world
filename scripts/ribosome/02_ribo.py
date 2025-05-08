#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches, animation

# ─── CONFIG ────────────────────────────────────────────────────────────────────
codon_width     = 0.5
num_codons      = 12
steps_per_codon = 15                   # sub-frames per codon (entry, dock, exit)
total_frames    = num_codons * steps_per_codon
interval_ms     = 100                  # ms between frames

# Example codon sequence and AA mapping:
codon_seqs = [
    "AUG", "GCU", "CGA", "UUU", "GAA", "CCU",
    "AGG", "UAC", "GGC", "AAA", "UUA", "GGG"
]
aa_map = {
    "AUG": "Met", "GCU": "Ala", "CGA": "Arg", "UUU": "Phe",
    "GAA": "Glu", "CCU": "Pro", "AGG": "Arg", "UAC": "Tyr",
    "GGC": "Gly", "AAA": "Lys", "UUA": "Leu", "GGG": "Gly"
}
# ───────────────────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(-1, 4)
ax.set_ylim(0, 2)
ax.axis("off")

# ─── Ribosome subunits ─────────────────────────────────────────────────────────
large_center = [1.0, 1.2]
small_center = [1.0, 0.5]
ax.add_patch(patches.Circle(large_center, 0.8, fill=True, alpha=0.3))
ax.add_patch(patches.Circle(small_center, 0.5, fill=True, alpha=0.3))

# ─── Site labels ───────────────────────────────────────────────────────────────
ax.annotate(
    "A-site",
    xy=(large_center[0], 0.35),          # near decoding site
    xytext=(large_center[0] + 0.6, 0.6),
    arrowprops=dict(arrowstyle="->"), fontsize=10
)
ax.annotate(
    "P-site",
    xy=(large_center[0] - 0.3, 1.15),     # inside large subunit
    xytext=(large_center[0] - 1.0, 1.5),
    arrowprops=dict(arrowstyle="->"), fontsize=10
)
ax.annotate(
    "E-site",
    xy=(large_center[0] + 1.0, 1.2),      # at exit tunnel
    xytext=(large_center[0] + 2.0, 1.5),
    arrowprops=dict(arrowstyle="->"), fontsize=10
)

# ─── mRNA + codon rectangles & labels ─────────────────────────────────────────
y_mRNA = 0.2
ax.plot([-1, 4], [y_mRNA, y_mRNA], lw=3)

# initial codon positions so first codon sits under A-site when aligned
target_center = large_center[0]
start_offset = target_center - codon_width/2
codon_positions = np.arange(num_codons) * codon_width + start_offset

codon_rects = []
codon_texts = []
for i, x in enumerate(codon_positions):
    rect = patches.Rectangle((x, y_mRNA - 0.05), codon_width, 0.1,
                             facecolor="lightgrey", edgecolor="black", alpha=0.5)
    ax.add_patch(rect)
    codon_rects.append(rect)
    txt = ax.text(
        x + codon_width/2, y_mRNA + 0.10, codon_seqs[i],
        ha="center", va="bottom", fontsize=10, fontweight="normal"
    )
    codon_texts.append(txt)

# ─── Peptide chain & AA labels ─────────────────────────────────────────────────
chain_scatter = ax.scatter([], [], s=100, color="green")
aa_texts      = []
for _ in range(num_codons):
    t = ax.text(0, 0, "", ha="center", va="center",
                fontsize=8, color="black", visible=False)
    aa_texts.append(t)

# ─── tRNA + AA + label ─────────────────────────────────────────────────────────
tRNA_scatter = ax.scatter([], [], marker="v", s=200, color="blue")
aa_scatter   = ax.scatter([], [], s=150, color="orange")
trna_label   = ax.text(0, 0, "tRNA", ha="center", va="bottom",
                       fontsize=12, color="blue", visible=False)

def init():
    empty = np.zeros((0,2))
    chain_scatter.set_offsets(empty)
    tRNA_scatter.set_offsets(empty)
    aa_scatter.set_offsets(empty)
    trna_label.set_visible(False)
    for txt in aa_texts:
        txt.set_visible(False)
    for rect, txt in zip(codon_rects, codon_texts):
        rect.set_facecolor("lightgrey")
        txt.set_fontweight("normal")
        rect.set_visible(True)
        txt.set_visible(True)
    return codon_rects + codon_texts + aa_texts + [chain_scatter, tRNA_scatter, aa_scatter, trna_label]

def animate(frame):
    empty = np.zeros((0,2))

    # 1) Slide mRNA (rectangles + labels)
    shift = (frame / steps_per_codon) * codon_width
    for i, (rect, txt) in enumerate(zip(codon_rects, codon_texts)):
        x = codon_positions[i] - shift
        rect.set_x(x)
        txt.set_x(x + codon_width/2)
        visible = (-1 < x < 4)
        rect.set_visible(visible)
        txt.set_visible(visible)
        rect.set_facecolor("lightgrey")
        txt.set_fontweight("normal")

    # 2) Highlight codon under A-site (true alignment)
    highlight_idx = None
    for i, pos in enumerate(codon_positions):
        x = pos - shift
        if x <= target_center <= x + codon_width:
            highlight_idx = i
            break
    if highlight_idx is not None:
        codon_rects[highlight_idx].set_facecolor("gold")
        codon_texts[highlight_idx].set_fontweight("bold")

    # 3) Grow peptide chain + AA labels
    codon_idx = frame // steps_per_codon
    aa_count  = min(codon_idx, num_codons)
    chain_x   = [large_center[0] + 0.8 + 0.1 + i*0.15 for i in range(aa_count)]
    chain_y   = [large_center[1]] * aa_count
    if aa_count:
        chain_scatter.set_offsets(np.column_stack([chain_x, chain_y]))
        for j in range(aa_count):
            aa_texts[j].set_text(aa_map[codon_seqs[j]])
            aa_texts[j].set_position((chain_x[j], chain_y[j]))
            aa_texts[j].set_visible(True)
    else:
        chain_scatter.set_offsets(empty)
        for txt in aa_texts:
            txt.set_visible(False)

    # 4) Animate tRNA entry/dock/exit
    sub = frame % steps_per_codon
    if codon_idx < num_codons:
        ent = steps_per_codon // 3
        dock = ent
        ex   = steps_per_codon - ent - dock
        if sub < ent:
            frac = sub / ent
            x = 4 - frac * (4 - large_center[0])
        elif sub < ent + dock:
            x = large_center[0]
        else:
            ef = (sub - ent - dock) / ex
            x = large_center[0] - ef * (large_center[0] + 1)
        y0 = y_mRNA + 0.15
        tRNA_scatter.set_offsets([[x, y0]])
        aa_scatter.set_offsets([[x, y0 + 0.25]])
        trna_label.set_position((x, y0 + 0.05))
        trna_label.set_visible(True)
    else:
        tRNA_scatter.set_offsets(empty)
        aa_scatter.set_offsets(empty)
        trna_label.set_visible(False)

    return codon_rects + codon_texts + aa_texts + [chain_scatter, tRNA_scatter, aa_scatter, trna_label]

anim = animation.FuncAnimation(
    fig, animate, init_func=init,
    frames=total_frames, interval=interval_ms, blit=True
)

plt.show()
