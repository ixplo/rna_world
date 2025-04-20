
import random
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib.animation import FuncAnimation

# Model parameters
NUCLEOTIDES = ['A', 'U', 'C', 'G']
RNA_LENGTH = 10
POPULATION_SIZE = 50
GENERATIONS = 30
MUTATION_RATE = 0.1
ENERGY_PER_RNA = 5
REPLICATION_COST = 2

def random_rna(length=RNA_LENGTH):
    return ''.join(random.choice(NUCLEOTIDES) for _ in range(length))

def mutate(rna, mutation_rate=MUTATION_RATE):
    return ''.join(
        n if random.random() > mutation_rate else random.choice(NUCLEOTIDES)
        for n in rna
    )

def is_ribozyme(rna):
    return "AUG" in rna

def replicate(rna, energy):
    if is_ribozyme(rna) and energy >= REPLICATION_COST:
        return mutate(rna), energy - REPLICATION_COST
    return None, energy

# Initialize population
population = [(random_rna(), ENERGY_PER_RNA) for _ in range(POPULATION_SIZE)]
history = []
sequence_history = []

for generation in range(GENERATIONS):
    next_gen = []
    current_sequences = []
    for rna, energy in population:
        current_sequences.append(rna)
        offspring, remaining_energy = replicate(rna, energy)
        if offspring:
            next_gen.append((rna, remaining_energy))
            next_gen.append((offspring, ENERGY_PER_RNA))
        else:
            if remaining_energy > 0:
                next_gen.append((rna, remaining_energy - 1))
    population = [pair for pair in next_gen if pair[1] > 0]
    history.append(len(population))
    sequence_history.append(Counter(current_sequences).most_common(5))

# Animation
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
fig.subplots_adjust(wspace=0.4)

line, = ax1.plot([], [], marker='o')
bar_container = ax2.bar([], [])

def init():
    ax1.set_xlim(0, GENERATIONS)
    ax1.set_ylim(0, max(history) + 10)
    ax1.set_title("Population Size Over Generations")
    ax1.set_xlabel("Generation")
    ax1.set_ylabel("Population")
    ax1.grid(True)
    return line, *bar_container

def update(frame):
    x = list(range(frame + 1))
    y = history[:frame + 1]
    line.set_data(x, y)

    ax2.clear()
    if sequence_history[frame]:
        labels, counts = zip(*sequence_history[frame])
        ax2.bar(labels, counts)
        ax2.set_title(f"Top 5 RNA (Generation {frame + 1})")
        ax2.set_ylabel("Frequency")
        ax2.set_xticks(range(len(labels)))
        ax2.set_xticklabels(labels, rotation=45, ha='right')
    return line, *bar_container

ani = FuncAnimation(fig, update, frames=GENERATIONS, init_func=init,
                    blit=False, interval=500, repeat=False)
plt.tight_layout()
plt.show()
