import random
import numpy as np
from scipy.ndimage import label
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, Button, RadioButtonGroup, Div, Span
from bokeh.plotting import figure
from bokeh.events import Tap
from collections import Counter
import copy

# ------- Simulation logic -------
GRID_SIZE = 20
INITIAL_MOLECULES = 200
GENERATIONS = 50
MAX_ENERGY = 5
REPLICATION_COST = 2
MUTATION_RATE = 0.1
MAX_TOTAL_MOLECULES = 5000
KERNEL_SIZE = 5

def random_rna(length=10):
    return ''.join(random.choice(['A','U','C','G']) for _ in range(length))

def is_ribozyme(rna):
    return 'AUG' in rna

def replicate(rna, energy):
    if is_ribozyme(rna) and energy >= REPLICATION_COST:
        offspring = ''.join(
            n if random.random() > MUTATION_RATE else random.choice(['A','U','C','G'])
            for n in rna
        )
        return offspring, energy - REPLICATION_COST
    return None, energy

def smooth_env(raw, ks):
    pad = ks // 2
    kernel = np.ones((ks, ks)) / (ks**2)
    padded = np.pad(raw, pad, mode='reflect')
    env = np.zeros_like(raw)
    for i in range(raw.shape[0]):
        for j in range(raw.shape[1]):
            env[i,j] = np.sum(padded[i:i+ks, j:j+ks] * kernel)
    return env

def run_simulation():
    N = GRID_SIZE
    raw_env = np.random.uniform(-1, 1, (N, N))
    ENV = smooth_env(raw_env, KERNEL_SIZE)
    clusters, _ = label(ENV > 0)
    grid = [[[] for _ in range(N)] for _ in range(N)]
    for _ in range(INITIAL_MOLECULES):
        x, y = random.randrange(N), random.randrange(N)
        grid[y][x].append((random_rna(), MAX_ENERGY))
    history = [copy.deepcopy(grid)]
    maps = []
    for _ in range(GENERATIONS):
        new_grid = [[[] for _ in range(N)] for _ in range(N)]
        for y in range(N):
            for x in range(N):
                for rna, energy in grid[y][x]:
                    energy = max(0, min(MAX_ENERGY, energy + ENV[y,x]))
                    offspring, rem = replicate(rna, energy)
                    if rem > 0:
                        new_grid[y][x].append((rna, rem))
                    if offspring and sum(len(c) for row in new_grid for c in row) < MAX_TOTAL_MOLECULES:
                        dx, dy = random.choice([(0,0),(1,0),(-1,0),(0,1),(0,-1)])
                        nx, ny = (x+dx)%N, (y+dy)%N
                        new_grid[ny][nx].append((offspring, MAX_ENERGY))
        grid = new_grid
        history.append(copy.deepcopy(grid))
        count_map = np.array([[len(cell) for cell in row] for row in grid])
        energy_map = np.array([[np.mean([e for (_,e) in cell]) if cell else 0 for cell in row] for row in grid])
        ribo_map = np.array([[sum(is_ribozyme(s) for (s,_) in cell)/len(cell) if cell else 0 for cell in row] for row in grid])
        maps.append({'count': count_map, 'energy': energy_map, 'ribo': ribo_map})
    return ENV, clusters, history, maps

# Run simulation
ENV, clusters, history, maps = run_simulation()
N = GRID_SIZE
gen_count = len(maps)

# Data sources
map_source = ColumnDataSource({'image': [maps[0]['count']]})
env_source = ColumnDataSource({'image': [ENV]})
metrics_source = ColumnDataSource(data=dict(
    gen=list(range(gen_count)),
    total=[int(np.sum(maps[g]['count'])) for g in range(gen_count)],
    ribo=[int(np.sum(maps[g]['ribo'] * maps[g]['count'])) for g in range(gen_count)],
    avgE=[float(np.mean(maps[g]['energy'][maps[g]['count']>0])) for g in range(gen_count)]
))

# Environment plot
env_fig = figure(title="Environment Levels", x_range=(0,N), y_range=(0,N), tools="wheel_zoom,reset", width=400, height=400)
env_fig.image('image', x=0, y=0, dw=N, dh=N, source=env_source, palette="RdBu11")

# Molecule count map
map_fig = figure(title="Count (Gen 0)", x_range=(0,N), y_range=(0,N), tools="tap,wheel_zoom,reset", width=400, height=400)
map_fig.image('image', x=0, y=0, dw=N, dh=N, source=map_source, palette="Viridis256")

# Detail panel
detail_div = Div(text="Click a cell for details", width=400, height=200)

def show_details(event):
    ix, iy = int(event.x), int(event.y)
    if 0 <= ix < N and 0 <= iy < N:
        gen = slider.value
        counts = [len(history[g][iy][ix]) for g in range(gen_count)]
        energies = [np.mean([e for (_,e) in history[g][iy][ix]]) if history[g][iy][ix] else 0 for g in range(gen_count)]
        ribo_fracs = [sum(is_ribozyme(s) for (s,_) in history[g][iy][ix]) / len(history[g][iy][ix]) if history[g][iy][ix] else 0 for g in range(gen_count)]
        seqs = [s for (s,_) in history[gen][iy][ix]]
        top = Counter(seqs).most_common(5)
        top_str = ", ".join(f"{s}×{c}" for s,c in top) or "None"
        detail_div.text = (
            f"<b>Cell ({ix},{iy}) at Gen {gen}</b><br>"
            f"Count: {counts[gen]}, AvgE: {energies[gen]:.2f}, Ribo%: {ribo_fracs[gen]*100:.1f}%<br>"
            f"<b>History Counts:</b> {counts}<br>"
            f"<b>History AvgE:</b> {[f'{e:.2f}' for e in energies]}<br>"
            f"<b>History Ribo%:</b> {[f'{r*100:.1f}' for r in ribo_fracs]}<br>"
            f"<b>Top seqs:</b> {top_str}"
        )

map_fig.on_event(Tap, show_details)

# Metrics figure
metrics_fig = figure(title="Metrics over Generations", width=800, height=300)
metrics_fig.line('gen','total', source=metrics_source, legend_label="Total")
metrics_fig.line('gen','ribo', source=metrics_source, legend_label="Ribo")
metrics_fig.line('gen','avgE', source=metrics_source, legend_label="AvgE")
vline = Span(location=0, dimension='height', line_dash='dashed', line_color='black')
metrics_fig.add_layout(vline)
metrics_fig.legend.location = "top_right"

# Controls
slider = Slider(start=0, end=gen_count-1, value=0, step=1, title="Generation", width=400)
map_buttons = RadioButtonGroup(labels=["count","energy","ribo"], active=0)

def update_map(attr, old, new):
    gen = slider.value
    key = map_buttons.labels[map_buttons.active]
    map_source.data = {'image': [maps[gen][key]]}
    map_fig.title.text = f"{key.capitalize()} (Gen {gen})"
    vline.location = gen

slider.on_change('value', update_map)
map_buttons.on_change('active', update_map)

# Play/Pause
callback_id = {'id': None}
def animate():
    gen = (slider.value + 1) % gen_count
    slider.value = gen

def play():
    if callback_id['id'] is None:
        callback_id['id'] = curdoc().add_periodic_callback(animate, 1000)

def pause():
    if callback_id['id']:
        curdoc().remove_periodic_callback(callback_id['id'])
        callback_id['id'] = None

play_button = Button(label="► Play")
pause_button = Button(label="❚❚ Pause")
play_button.on_click(play)
pause_button.on_click(pause)

# Layout
layout = column(
    row(env_fig, map_fig, detail_div),
    row(slider, play_button, pause_button, map_buttons),
    metrics_fig
)

curdoc().add_root(layout)
curdoc().title = "RNA Spatial Simulation Dashboard v2 (Fixed)"