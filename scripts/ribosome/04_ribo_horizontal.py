import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches, animation

speed = 1.0
belt_dx = 0.05 * speed
fall_vy = 0.025 * speed
interval = int(100 / speed)

num_codons = 20
codon_w, codon_h = 1.0, 0.3
spacing = codon_w
teeth = 3
colors = plt.cm.tab20(np.linspace(0,1,num_codons))
belt_len = num_codons * spacing
y_start = 2.0
belt_y = codon_h / 2
spawn_interval = int(spacing / belt_dx)

offset = 0.15
spawn_x = belt_len + (y_start - belt_y) * (belt_dx / fall_vy) - codon_w / 2 + offset
landing_x = belt_len / 2 - codon_w / 2 + offset
N_land = int((y_start - belt_y) / fall_vy)
vxfly = (spawn_x - landing_x) / N_land

fig, ax = plt.subplots(figsize=(10, 3))
ax.set_aspect('equal')
ax.axis('off')
ax.set_xlim(0, belt_len)
ax.set_ylim(-1, y_start + 1)

codons = []
for i in range(num_codons):
    x0 = belt_len - i * spacing
    verts = [(x0, -codon_h/2), (x0, codon_h/2)]
    for j in range(teeth):
        xm = x0 + (j + 0.5) * codon_w / teeth
        xb = x0 + (j + 1) * codon_w / teeth
        verts += [(xm, codon_h/4), (xb, codon_h/2)]
    verts.append((x0 + codon_w, -codon_h/2))
    poly = patches.Polygon(verts, closed=True, facecolor=colors[i],
                           edgecolor='black', zorder=1)
    ax.add_patch(poly)
    codons.append(poly)

shapes = []

def get_color_for_position(x):
    for codon in codons:
        codon_xs = codon.get_xy()[:, 0]
        if np.min(codon_xs) <= x <= np.max(codon_xs):
            return codon.get_facecolor()
    return 'gray'

def rotate(verts, angle_deg, center):
    angle = np.radians(angle_deg)
    cx, cy = center
    rot = []
    for x, y in verts:
        dx, dy = x - cx, y - cy
        x_new = cx + dx * np.cos(angle) - dy * np.sin(angle)
        y_new = cy + dx * np.sin(angle) + dy * np.cos(angle)
        rot.append((x_new, y_new))
    return rot

def make_trna_patch(xc, yc, color='gray', angle=0):
    w = codon_w
    verts = []
    left = xc - w/2
    verts.append((left, yc))
    for j in range(teeth):
        x1 = left + j * (w/teeth)
        xm = x1 + (w/teeth)/2
        x2 = x1 + w/teeth
        verts += [(xm, yc - codon_h/4), (x2, yc)]
    verts.append((xc + w/2, yc))

    stem_h = codon_h * 1.2
    stem_w = w * 0.4
    bar_h = codon_h * 0.4
    sw = stem_w / 2
    verts += [
        (xc + w/2, yc + stem_h),
        (xc + sw, yc + stem_h),
        (xc + sw, yc + stem_h + bar_h),
        (xc - sw, yc + stem_h + bar_h),
        (xc - sw, yc + stem_h),
        (xc - w/2, yc + stem_h),
        (xc - w/2, yc)
    ]
    if angle != 0:
        verts = rotate(verts, angle, center=(xc, yc))
    return patches.Polygon(verts, closed=True, facecolor=color,
                           edgecolor='black', zorder=5)

def spawn_shape():
    predicted_shift = belt_dx * N_land
    target_x = landing_x + predicted_shift + codon_w / 2
    target_color = get_color_for_position(target_x)
    p = make_trna_patch(spawn_x, y_start, target_color)
    ax.add_patch(p)
    shapes.append({
        'patch': p,
        'state': 'flying',
        'x': spawn_x,
        'y': y_start,
        'color': target_color,
        'angle': 0
    })

def update(frame):
    if frame and frame % spawn_interval == 0:
        spawn_shape()

    for poly in codons:
        verts = poly.get_xy()
        verts[:, 0] -= belt_dx
        if verts[0, 0] < -codon_w:
            verts[:, 0] += belt_len
        poly.set_xy(verts)

    for s in shapes:
        if s['state'] == 'flying':
            s['x'] -= vxfly
            s['y'] -= fall_vy
            if s['y'] <= belt_y:
                s['state'] = 'landed'
                s['y'] = belt_y
                s['x'] = landing_x
        elif s['state'] == 'landed':
            s['x'] -= belt_dx
        elif s['state'] == 'detached':
            s['y'] -= fall_vy
            s['angle'] += 2  # rotate counter-clockwise

        s['patch'].set_xy(make_trna_patch(s['x'], s['y'], s['color'], s['angle']).get_xy())

    landed = [s for s in shapes if s['state'] == 'landed']
    if len(landed) > 2:
        landed[0]['state'] = 'detached'

    shapes[:] = [s for s in shapes if not (s['state'] == 'detached' and s['y'] < -1)]

    return codons + [s['patch'] for s in shapes]

anim = animation.FuncAnimation(fig, update, frames=1000,
                               interval=interval, blit=True)
plt.show()