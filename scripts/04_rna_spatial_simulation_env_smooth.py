import random
import matplotlib.pyplot as plt
import numpy as np

# Spatial RNA-world simulation with smoothed environment visualization

GRID_SIZE = 20
INITIAL_MOLECULES = 200
MAX_ENERGY = 5
REPLICATION_COST = 2
MUTATION_RATE = 0.1
GENERATIONS = 50
MAX_TOTAL_MOLECULES = 5000

# Colormaps
CMAP_COUNT = 'viridis'
CMAP_ENV = 'coolwarm'

NUCLEOTIDES = ['A', 'U', 'C', 'G']
RNA_LENGTH = 10

def random_rna(length=RNA_LENGTH):
    return ''.join(random.choice(NUCLEOTIDES) for _ in range(length))

def mutate(rna, mutation_rate=MUTATION_RATE):
    return ''.join(n if random.random() > mutation_rate else random.choice(NUCLEOTIDES) for n in rna)

def is_ribozyme(rna):
    return "AUG" in rna

def replicate(rna, energy):
    if is_ribozyme(rna) and energy >= REPLICATION_COST:
        return mutate(rna), energy - REPLICATION_COST
    return None, energy

# Generate raw environment grid and apply smoothing for contiguous regions
raw_env = np.random.uniform(-1, 1, (GRID_SIZE, GRID_SIZE))
kernel_size = 5
kernel = np.ones((kernel_size, kernel_size)) / (kernel_size**2)
pad = kernel_size // 2
padded = np.pad(raw_env, pad, mode='reflect')
ENV_GRID = np.zeros_like(raw_env)
for i in range(GRID_SIZE):
    for j in range(GRID_SIZE):
        ENV_GRID[i, j] = np.sum(padded[i:i+kernel_size, j:j+kernel_size] * kernel)

# Initialize molecule grid
grid = [[[] for _ in range(GRID_SIZE)] for _ in range(GRID_SIZE)]
for _ in range(INITIAL_MOLECULES):
    x, y = random.randrange(GRID_SIZE), random.randrange(GRID_SIZE)
    grid[y][x].append((random_rna(), MAX_ENERGY))

# Set up interactive plot: two side-by-side maps
plt.ion()
fig, (env_ax, mol_ax) = plt.subplots(1, 2, figsize=(12, 6))

env_im = env_ax.imshow(ENV_GRID, cmap=CMAP_ENV, interpolation='nearest')
env_ax.set_title('Smoothed Environment')
fig.colorbar(env_im, ax=env_ax, fraction=0.046, pad=0.04, label='Env Value')

counts = np.array([[len(grid[y][x]) for x in range(GRID_SIZE)] for y in range(GRID_SIZE)])
mol_im = mol_ax.imshow(counts, cmap=CMAP_COUNT, interpolation='nearest')
mol_ax.set_title('Molecule Density')
fig.colorbar(mol_im, ax=mol_ax, fraction=0.046, pad=0.04, label='Count')

stats_text = fig.text(0.5, 0.02, "", ha='center', va='bottom', fontsize=10)

for gen in range(GENERATIONS):
    offspring_count = 0
    new_grid = [[[] for _ in range(GRID_SIZE)] for _ in range(GRID_SIZE)]
    for y in range(GRID_SIZE):
        for x in range(GRID_SIZE):
            for rna, energy in grid[y][x]:
                # Adjust energy by smoothed environment
                energy = max(0, min(MAX_ENERGY, energy + ENV_GRID[y][x]))
                offspring, rem_energy = replicate(rna, energy)
                if rem_energy > 0:
                    new_grid[y][x].append((rna, rem_energy))
                if offspring and sum(len(cell) for row in new_grid for cell in row) < MAX_TOTAL_MOLECULES:
                    dx, dy = random.choice([(0,0), (1,0), (-1,0), (0,1), (0,-1)])
                    nx, ny = (x + dx) % GRID_SIZE, (y + dy) % GRID_SIZE
                    new_grid[ny][nx].append((offspring, MAX_ENERGY))
                    offspring_count += 1
    grid = new_grid

    # Compute statistics
    molecules = [(rna, energy) for row in grid for cell in row for (rna, energy) in cell]
    total = len(molecules)
    ribo_count = sum(1 for rna, _ in molecules if is_ribozyme(rna))
    avg_energy = np.mean([energy for _, energy in molecules]) if total > 0 else 0
    unique_seq = len(set(rna for rna, _ in molecules))
    avg_env = (ENV_GRID * np.array([[len(grid[y][x]) for x in range(GRID_SIZE)] 
             for y in range(GRID_SIZE)])).sum() / total if total > 0 else 0

    # Update visualizations
    counts = np.array([[len(grid[y][x]) for x in range(GRID_SIZE)] for y in range(GRID_SIZE)])
    mol_im.set_data(counts)
    mol_ax.set_title(f'Generation {gen}')
    stats = (
        f"Gen: {gen} | Total: {total} | Ribozymes: {ribo_count} | "
        f"Avg Energy: {avg_energy:.2f} | Unique Seq: {unique_seq} | "
        f"Offspring: {offspring_count} | Avg Env: {avg_env:.2f}"
    )
    stats_text.set_text(stats)
    fig.canvas.draw()
    plt.pause(1.5)

plt.ioff()
plt.show()