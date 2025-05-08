#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches, animation

# ─── CONFIG ────────────────────────────────────────────────────────────────────
codon_width       = 0.5
num_codons        = 12
steps_per_codon   = 15                   # sub-frames per codon
total_frames      = num_codons * steps_per_codon
interval_ms       = 100                  # ms between frames

# You can supply your real codon sequences here:
codon_seqs = [
    "AUG", "GCU", "CGA", "UUU", "GAA", "CCU",
    "AGG", "UAC", "GGC", "AAA", "UUA", "GGG"
]
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

# ─── mRNA line + codons ────────────────────────────────────────────────────────
y_mRNA = 0.2
ax.plot([-1, 4], [y_mRNA, y_mRNA], lw=3)

codon_positions = np.arange(num_codons) * codon_width
codon_rects = []
codon_texts = []
for i, x in enumerate(codon_positions):
    # rectangle
    rect = patches.Rectangle((x, y_mRNA - 0.05), codon_width, 0.1, alpha=0.5)
    ax.add_patch(rect)
    codon_rects.append(rect)
    # text label
    txt = ax.text(
        x + codon_width/2,      # center of rectangle
        y_mRNA + 0.10,          # just above the mRNA line
        codon_seqs[i],          # sequence
        ha="center", va="bottom",
        fontsize=10, color="black"
    )
    codon_texts.append(txt)

# ─── Growing peptide chain ─────────────────────────────────────────────────────
chain_scatter = ax.scatter([], [], s=100, color="green")

# ─── tRNA + amino acid + label ─────────────────────────────────────────────────
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
    for rect, txt in zip(codon_rects, codon_texts):
        rect.set_visible(True)
        txt.set_visible(True)
    return codon_rects + codon_texts + [chain_scatter, tRNA_scatter, aa_scatter, trna_label]

def animate(frame):
    empty = np.zeros((0,2))

    # 1) Slide the mRNA (rectangles + labels) left
    shift = (frame / steps_per_codon) * codon_width
    for i, (rect, txt) in enumerate(zip(codon_rects, codon_texts)):
        x = codon_positions[i] - shift
        rect.set_x(x)
        txt.set_x(x + codon_width/2)
        visible = (-1 < x < 4)
        rect.set_visible(visible)
        txt.set_visible(visible)

    # 2) Grow the peptide chain
    codon_idx = frame // steps_per_codon
    aa_count  = min(codon_idx, num_codons)
    chain_x   = [large_center[0] + 0.8 + 0.1 + i*0.15 for i in range(aa_count)]
    chain_y   = [large_center[1]] * aa_count
    if aa_count:
        chain_scatter.set_offsets(np.column_stack([chain_x, chain_y]))
    else:
        chain_scatter.set_offsets(empty)

    # 3) Animate tRNA entry, dock, & exit
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

    return codon_rects + codon_texts + [chain_scatter, tRNA_scatter, aa_scatter, trna_label]

anim = animation.FuncAnimation(
    fig, animate, init_func=init,
    frames=total_frames, interval=interval_ms, blit=True
)

plt.show()
