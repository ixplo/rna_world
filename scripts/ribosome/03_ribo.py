#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches, animation

# ─── CONFIG ──────────────────────────────────────────────────────────────────
codon_width     = 0.5
num_codons      = 12
steps_per_codon = 15                   # frames per codon
total_frames    = num_codons * steps_per_codon
interval_ms     = 100                  # ms between frames
exit_angle_deg  = 45                   # peptide chain angle
spacing         = 0.15                 # distance between AA on chain

# Codon sequences and AA mapping
codon_seqs = [
    "AUG", "GCU", "CGA", "UUU", "GAA", "CCU",
    "AGG", "UAC", "GGC", "AAA", "UUA", "GGG"
]
aa_map = {
    "AUG": "Met", "GCU": "Ala", "CGA": "Arg", "UUU": "Phe",
    "GAA": "Glu", "CCU": "Pro", "AGG": "Arg", "UAC": "Tyr",
    "GGC": "Gly", "AAA": "Lys", "UUA": "Leu", "GGG": "Gly"
}
# Colors for chain spheres
colors = plt.cm.tab20(np.linspace(0, 1, num_codons))

# Pre-calc
exit_angle = np.deg2rad(exit_angle_deg)

# Set up figure
fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(-1, 4)
ax.set_ylim(0, 2)
ax.axis('off')

# Ribosome subunits
large_center = np.array([1.0, 1.2])
large_radius = 0.8
small_center = np.array([1.0, 0.5])
ax.add_patch(patches.Circle(large_center, large_radius, fill=True, alpha=0.3))
ax.add_patch(patches.Circle(small_center, 0.5, fill=True, alpha=0.3))

# mRNA baseline
y_mRNA = 0.2
ax.plot([-1, 4], [y_mRNA, y_mRNA], lw=3)

# Label A, P, E site positions relative to ribosome subunits
A_site_x = small_center[0]
A_site_y = small_center[1] + 0.05
P_site_x = large_center[0] - 0.2
P_site_y = large_center[1] - 0.1
E_site_x = large_center[0] + large_radius * np.cos(exit_angle)
E_site_y = large_center[1] + large_radius * np.sin(exit_angle)


# mRNA + codons
target_third = 2  # index of codon under A-site at start
start_offset = large_center[0] - codon_width/2 - target_third * codon_width
codon_positions = np.arange(num_codons) * codon_width + start_offset
codon_rects, codon_texts = [], []
for i, x in enumerate(codon_positions):
    rect = patches.Rectangle((x, y_mRNA - 0.05), codon_width, 0.1,
                             facecolor='lightgrey', edgecolor='black', alpha=0.5)
    ax.add_patch(rect)
    codon_rects.append(rect)
    txt = ax.text(x + codon_width/2, y_mRNA + 0.1, codon_seqs[i],
                  ha='center', va='bottom', fontsize=10)
    codon_texts.append(txt)

#
# Peptide chain: circles and texts
chain_scatter = ax.scatter([], [], s=100)  # will not be used for drawing anymore, but kept for return
aa_texts = [ax.text(0, 0, '', ha='center', va='center', fontsize=8,
                    color='black', visible=False)
            for _ in range(num_codons)]
# List to hold chain circles (patches)
chain_circles = []

# tRNA conveyors and attached AA
tRNA_scats = [ax.scatter([], [], marker='v', s=150, color='blue') for _ in range(3)]
aa_scats   = [ax.scatter([], [], s=100, color='orange')   for _ in range(3)]

# INITIALIZE
def init():
    empty = np.zeros((0, 2))
    chain_scatter.set_offsets(empty)
    chain_scatter.set_color([])
    for scat in tRNA_scats + aa_scats:
        scat.set_offsets(empty)
    for txt in aa_texts:
        txt.set_visible(False)
    for rect, txt in zip(codon_rects, codon_texts):
        rect.set_facecolor('lightgrey')
        txt.set_fontweight('normal')
    # Remove all chain circles if any exist
    global chain_circles
    for circ in chain_circles:
        circ.remove()
    chain_circles = []
    return codon_rects + codon_texts + aa_texts + tRNA_scats + aa_scats

# ANIMATE
def animate(frame):
    empty = np.zeros((0, 2))
    shift = (frame / steps_per_codon) * codon_width

    # Slide mRNA and reset highlights
    for rect, txt, pos in zip(codon_rects, codon_texts, codon_positions):
        x = pos - shift
        rect.set_x(x)
        txt.set_x(x + codon_width/2)
        vis = -1 < x < 4
        rect.set_visible(vis)
        txt.set_visible(vis)
        rect.set_facecolor('lightgrey')
        txt.set_fontweight('normal')

    # Highlight current codon under A-site
    a_x = A_site_x
    idx = next((i for i, pos in enumerate(codon_positions)
                if pos - shift <= a_x <= pos - shift + codon_width), None)
    if idx is not None:
        codon_rects[idx].set_facecolor('gold')
        codon_texts[idx].set_fontweight('bold')

    # Determine current codon index for peptide chain
    codon_idx = frame // steps_per_codon
    aa_count = min(codon_idx, num_codons)

    # Remove old chain circles
    global chain_circles
    for circ in chain_circles:
        circ.remove()
    chain_circles = []

    chain_pos = []
    for i in range(aa_count):
        rad = spacing * (aa_count - 1 - i)
        x = A_site_x
        y = A_site_y + 0.3 + rad
        chain_pos.append((x, y))
    if aa_count:
        # Draw circles and set texts
        for j, (x, y) in enumerate(chain_pos):
            circ = patches.Circle((x, y), 0.08, facecolor=colors[j], edgecolor='black', zorder=5)
            ax.add_patch(circ)
            chain_circles.append(circ)
            aa_texts[j].set_text(aa_map[codon_seqs[j]])
            aa_texts[j].set_position((x, y))
            aa_texts[j].set_visible(True)
        # Hide unused texts
        for k in range(aa_count, num_codons):
            aa_texts[k].set_visible(False)
    else:
        for txt in aa_texts:
            txt.set_visible(False)

    # chain_scatter is not used for drawing, but keep it empty for blitting
    chain_scatter.set_offsets(empty)
    chain_scatter.set_color([])

    # tRNA conveyors: A-site = current, P-site = current-1, E-site = current-2
    offsets = [0, -1, -2]
    sites = [(A_site_x, A_site_y), (P_site_x, P_site_y), (E_site_x, E_site_y)]
    for k, off in enumerate(offsets):
        t_idx = frame // steps_per_codon + off
        if 0 <= t_idx < num_codons:
            progress = (frame % steps_per_codon) / steps_per_codon
            if k == 0:  # incoming
                x_from, y_from = A_site_x + 2.0, A_site_y + 1.0
                x_to, y_to = A_site_x, A_site_y
                x_c = x_from + (x_to - x_from) * progress
                y_c = y_from + (y_to - y_from) * progress

                # Set the color of incoming amino acid
                if t_idx < num_codons:
                    aa_scats[k].set_color(colors[t_idx])
            elif k == 1:  # current attached, align to highlighted codon
                if idx is not None:
                    codon_x = codon_positions[idx] - shift
                    x_c = codon_x + codon_width / 2
                    y_c = A_site_y
                else:
                    x_c, y_c = A_site_x, A_site_y
            else:  # exiting
                x_from, y_from = P_site_x, P_site_y
                x_to, y_to = P_site_x - 2.0, P_site_y + 1.0
                x_c = x_from + (x_to - x_from) * progress
                y_c = y_from + (y_to - y_from) * progress
            tRNA_scats[k].set_offsets([[x_c, y_c]])
            aa_scats[k].set_offsets([[x_c, y_c + 0.2]])
        else:
            tRNA_scats[k].set_offsets(empty)
            aa_scats[k].set_offsets(empty)

    return codon_rects + codon_texts + aa_texts + tRNA_scats + aa_scats + chain_circles

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=total_frames, interval=interval_ms, blit=True)
plt.show()