import argparse
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.animation import FuncAnimation
from collections import Counter
import mplcursors
from scipy.ndimage import label

def parse_args():
    parser = argparse.ArgumentParser(description='Spatial RNA-world simulation with interactive visualization.')
    parser.add_argument('--grid_size', type=int, default=20, help='Grid dimension (NxN)')
    parser.add_argument('--initial', type=int, default=200, help='Initial number of RNA molecules')
    parser.add_argument('--generations', type=int, default=50, help='Number of generations')
    parser.add_argument('--max_energy', type=int, default=5)
    parser.add_argument('--rep_cost', type=int, default=2)
    parser.add_argument('--mutation', type=float, default=0.1)
    parser.add_argument('--max_molecules', type=int, default=5000)
    parser.add_argument('--kernel', type=int, default=5, help='Smoothing kernel size for environment')
    parser.add_argument('--export', action='store_true', help='Export animation to mp4')
    return parser.parse_args()

def random_rna(length=10):
    return ''.join(random.choice(['A','U','C','G']) for _ in range(length))

def is_ribozyme(seq):
    return 'AUG' in seq

def replicate(rna, energy, rep_cost, mutation):
    if is_ribozyme(rna) and energy >= rep_cost:
        offspring = ''.join(n if random.random()>mutation else random.choice(['A','U','C','G']) for n in rna)
        return offspring, energy - rep_cost
    return None, energy

def smooth_env(raw, kernel_size):
    pad = kernel_size//2
    kernel = np.ones((kernel_size, kernel_size)) / (kernel_size**2)
    padded = np.pad(raw, pad, mode='reflect')
    out = np.zeros_like(raw)
    for i in range(raw.shape[0]):
        for j in range(raw.shape[1]):
            out[i,j] = np.sum(padded[i:i+kernel_size, j:j+kernel_size] * kernel)
    return out

def run_simulation(args):
    N = args.grid_size
    raw_env = np.random.uniform(-1, 1, (N, N))
    ENV = smooth_env(raw_env, args.kernel)
    env_mask = ENV > 0
    clusters, num_clusters = label(env_mask)
    sizes = [np.sum(clusters==i) for i in range(1, num_clusters+1)]
    avg_cluster = np.mean(sizes) if sizes else 0
    grid = [[[] for _ in range(N)] for _ in range(N)]
    for _ in range(args.initial):
        x,y = random.randrange(N), random.randrange(N)
        grid[y][x].append((random_rna(), args.max_energy))
    maps = []
    metrics = {'total':[], 'ribo':[], 'avgE':[], 'unique':[], 'avgEnv':[]}
    for gen in range(args.generations):
        new_grid = [[[] for _ in range(N)] for _ in range(N)]
        for i in range(N):
            for j in range(N):
                for seq, energy in grid[i][j]:
                    energy = max(0, min(args.max_energy, energy + ENV[i,j]))
                    offspring, rem = replicate(seq, energy, args.rep_cost, args.mutation)
                    if rem > 0:
                        new_grid[i][j].append((seq, rem))
                    if offspring and sum(len(c) for row in new_grid for c in row) < args.max_molecules:
                        dx,dy = random.choice([(0,0),(1,0),(-1,0),(0,1),(0,-1)])
                        ni, nj = (i+dy)%N, (j+dx)%N
                        new_grid[ni][nj].append((offspring, args.max_energy))
        grid = new_grid
        count_map = np.array([[len(cell) for cell in row] for row in grid])
        energ_map = np.array([[np.mean([e for (_,e) in cell]) if cell else 0 for cell in row] for row in grid])
        ribo_map = np.array([[sum(is_ribozyme(s) for (s,_) in cell)/len(cell) if cell else 0 for cell in row] for row in grid])
        maps.append({'count':count_map, 'energy':energ_map, 'ribo':ribo_map})
        all_cells = [(s,e) for row in grid for cell in row for (s,e) in cell]
        totals = len(all_cells)
        metrics['total'].append(totals)
        metrics['ribo'].append(sum(is_ribozyme(s) for s,_ in all_cells))
        metrics['avgE'].append(np.mean([e for _,e in all_cells]) if all_cells else 0)
        metrics['unique'].append(len(set(s for s,_ in all_cells)))
        metrics['avgEnv'].append((ENV*count_map).sum()/totals if totals else 0)
    return ENV, clusters, num_clusters, avg_cluster, maps, metrics

def visualize(args, ENV, clusters, num_clusters, avg_cluster, maps, metrics):
    N = args.grid_size
    fig = plt.figure(constrained_layout=True, figsize=(10,8))
    gs = fig.add_gridspec(3, 3)
    ax_env = fig.add_subplot(gs[0, 0])
    ax_map = fig.add_subplot(gs[0, 1:])
    ax_metrics = fig.add_subplot(gs[1:, 0:2])
    ax_slider = fig.add_subplot(gs[2,2])

    im_env = ax_env.imshow(ENV, cmap='coolwarm', interpolation='nearest')
    ax_env.set_title('Environment')
    plt.colorbar(im_env, ax=ax_env)

    im_map = ax_map.imshow(maps[0]['count'], cmap='viridis', interpolation='nearest')
    ax_map.set_title('Count (Gen 0)')
    plt.colorbar(im_map, ax=ax_map)

    gens = list(range(len(metrics['total'])))
    ax_metrics.plot(gens, metrics['total'], label='Total')
    ax_metrics.plot(gens, metrics['ribo'],  label='Ribo')
    ax_metrics.plot(gens, metrics['avgE'],  label='AvgE')
    l_line = ax_metrics.axvline(0, color='k', linestyle='--')
    ax_metrics.legend(loc='upper right')
    ax_metrics.set_xlabel('Generation')

    slider = Slider(ax_slider, 'Gen', 0, len(maps)-1, valinit=0, valstep=1)
    ax_radio = fig.add_axes([0.8, 0.5, 0.1, 0.15])
    radio = RadioButtons(ax_radio, ('count','energy','ribo'))
    current = 'count'

    def update(gen):
        im_map.set_data(maps[gen][current])
        ax_map.set_title(f'{current} (Gen {gen})')
        # Prevent recursion by disabling slider events
        slider.eventson = False
        slider.set_val(gen)
        slider.eventson = True
        l_line.set_xdata([gen, gen])
        fig.canvas.draw_idle()

    slider.on_changed(lambda val: update(int(val)))
    radio.on_clicked(lambda label: (model := label, update(int(slider.val))))

    ax_play = fig.add_axes([0.8, 0.7, 0.1, 0.05])
    btn_play = Button(ax_play, 'Play/Pause')
    playing = {'status': False}
    anim = {'obj': None}

    def animate(frame):
        update(frame)

    def on_play(event):
        playing['status'] = not playing['status']
        if playing['status']:
            anim['obj'] = FuncAnimation(fig, animate, frames=len(maps), interval=1000, repeat=False)
            if args.export:
                anim['obj'].save('simulation.mp4', writer='ffmpeg')

    btn_play.on_clicked(on_play)

    cursor = mplcursors.cursor(im_map, hover=True)
    @cursor.connect("add")
    def on_add(sel):
        row, col = sel.index
        cnt = maps[int(slider.val)]['count'][row, col]
        avgE = maps[int(slider.val)]['energy'][row, col]
        ribo_frac = maps[int(slider.val)]['ribo'][row, col]
        cl_id = clusters[row, col]
        info = (
            f"Cell ({row},{col})\n"
            f"Count: {cnt}\n"
            f"Avg Energy: {avgE:.2f}\n"
            f"Ribo %: {ribo_frac*100:.1f}\n"
            f"Env Cluster: {cl_id}"
        )
        sel.annotation.set_text(info)

    plt.show()

def main():
    args = parse_args()
    ENV, clusters, num_clusters, avg_cluster, maps, metrics = run_simulation(args)
    print(f'Env clusters: {num_clusters}, avg size={avg_cluster:.2f}')
    visualize(args, ENV, clusters, num_clusters, avg_cluster, maps, metrics)

if __name__ == '__main__':
    main()