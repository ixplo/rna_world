import random
import matplotlib.pyplot as plt
import numpy as np

# Spatial RNA-world simulation with per-generation animation and detailed statistics
# Updated colormap: viridis

GRID_SIZE = 20
INITIAL_MOLECULES = 200
MAX_ENERGY = 5
REPLICATION_COST = 2
MUTATION_RATE = 0.1
GENERATIONS = 50
MAX_TOTAL_MOLECULES = 5000
CMAP = 'viridis'  # perceptually uniform colormap, avoids white-on-white

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

# Initialize grid
grid = [[[] for _ in range(GRID_SIZE)] for _ in range(GRID_SIZE)]
for _ in range(INITIAL_MOLECULES):
    x, y = random.randrange(GRID_SIZE), random.randrange(GRID_SIZE)
    grid[y][x].append((random_rna(), MAX_ENERGY))

# Set up plot
plt.ion()
fig, ax = plt.subplots(figsize=(6, 6))

# Initial display
counts = np.array([[len(grid[y][x]) for x in range(GRID_SIZE)] for y in range(GRID_SIZE)])
im = ax.imshow(counts, cmap=CMAP, interpolation='nearest')
cbar = fig.colorbar(im, ax=ax, label='Molecule count')

# Text box for statistics
stats_text = fig.text(0.02, 0.95, "", va='top', fontsize=10)

for gen in range(GENERATIONS):
    offspring_count = 0
    new_grid = [[[] for _ in range(GRID_SIZE)] for _ in range(GRID_SIZE)]
    for y in range(GRID_SIZE):
        for x in range(GRID_SIZE):
            for rna, energy in grid[y][x]:
                offspring, rem_energy = replicate(rna, energy)
                if rem_energy > 0:
                    new_grid[y][x].append((rna, rem_energy))
                if offspring and sum(len(cell) for row in new_grid for cell in row) < MAX_TOTAL_MOLECULES:
                    dx, dy = random.choice([(0,0), (1,0), (-1,0), (0,1), (0,-1)])
                    nx, ny = (x + dx) % GRID_SIZE, (y + dy) % GRID_SIZE
                    new_grid[ny][nx].append((offspring, MAX_ENERGY))
                    offspring_count += 1
    grid = new_grid

    # Compute detailed statistics
    molecules = [(rna, energy) for row in grid for cell in row for (rna, energy) in cell]
    total = len(molecules)
    ribo_count = sum(1 for rna, _ in molecules if is_ribozyme(rna))
    avg_energy = np.mean([energy for _, energy in molecules]) if total > 0 else 0
    unique_seq = len(set(rna for rna, _ in molecules))

    # Update plot
    counts = np.array([[len(grid[y][x]) for x in range(GRID_SIZE)] for y in range(GRID_SIZE)])
    im.set_data(counts)
    ax.set_title(f'Generation {gen}')
    stats = (
        f"Total molecules: {total}\n"
        f"Ribozymes: {ribo_count}\n"
        f"Avg energy: {avg_energy:.2f}\n"
        f"Unique sequences: {unique_seq}\n"
        f"Offspring created: {offspring_count}"
    )
    stats_text.set_text(stats)
    fig.canvas.draw()
    plt.pause(1.5)

plt.ioff()
plt.show()