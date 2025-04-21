# RNA World Simulation

This repository contains simple Python simulations that demonstrate how life could have emerged from the RNA world — an early stage of evolution before DNA and proteins.

We simulate:
- Self-replicating RNA molecules (ribozymes)
- Evolution over generations
- Mutations and selection
- Spatial behavior in a 2D environment

## Demo

Welcome to this RNA World simulation showcase.

🧬 We begin by running the first script:
01_rna_world_simulation.py.
This simulates a population of RNA sequences that replicate, mutate, and evolve over time in a one-dimensional space.
Watch how the population changes with each generation — and how some sequences begin to dominate.

🌐 Next, we move to 02_rna_spatial_simulation_colormap.py.
This introduces a two-dimensional environment.
Each cell on the grid can contain RNA, and replication now depends on local surroundings.
The color map shows the density of sequences across the grid.

🌱 In 03_rna_spatial_simulation_with_env_viz.py, we introduce environmental factors.
Some areas are richer in energy, and RNA thrives better there.
The overlay shows both RNA presence and the underlying environment.

🌄 Now, 04_rna_spatial_simulation_env_smooth.py smooths the environmental transitions.
This creates a more natural terrain — favoring gradual adaptations instead of abrupt changes.

🕹️ Script 05_rna_spatial_simulation_interactive.py brings interactivity.
We can now adjust parameters like mutation rate or energy per cell and observe the effects live.

📊 06_rna_spatial_simulation_with_statistics.py collects useful metrics.
We track population size, diversity, and evolutionary patterns.
These graphs help us understand how life-like behaviors emerge over time.

📈 Finally, 07_rna_spatial_simulation_dashboard_bokeh.py launches a full dashboard.
Built with Bokeh, it lets us explore regions, zoom in on local populations, and visualize history interactively.

Watch the YouTube video: [https://youtu.be/OhSRjQT50h4]

## Files

- `scripts/00_matplotlib_warmup.py`: Warm-up script to trigger matplotlib font cache
- `scripts/01_rna_world_simulation.py`: Basic evolutionary model (1D population)
- `scripts/02_rna_spatial_simulation_colormap.py`: 2D simulation with color map visualization
- `scripts/03_rna_spatial_simulation_with_env_viz.py`: Adds environmental visualization
- `scripts/04_rna_spatial_simulation_env_smooth.py`: Smoothing for environmental variables
- `scripts/05_rna_spatial_simulation_interactive.py`: Interactive visualization
- `scripts/06_rna_spatial_simulation_with_statistics.py`: Collects and visualizes statistics
- `scripts/07_rna_spatial_simulation_dashboard_bokeh.py`: Bokeh dashboard for interactive analysis
- `rna_world_youtube_script.txt`: Full video script
- `rna_world_subtitles.srt`: Sample subtitles
- `rna_world_youtube_description.txt`: Optimized YouTube description

## Requirements

- Python 3.9+
- `matplotlib`
- `numpy`

Install with:

```bash
pip install matplotlib numpy
```

Using venv:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib scipy mplcursors ffmpeg-python
pip install bokeh
```

## How to run

### Run 00_matplotlib_warmup.py

```bash
python3 scripts/00_matplotlib_warmup.py
```

### Run 01_rna_world_simulation.py

```bash
python3 scripts/01_rna_world_simulation.py
```

### Run 02_rna_spatial_simulation_colormap.py

```bash
python3 scripts/02_rna_spatial_simulation_colormap.py
```

### Run 03_rna_spatial_simulation_with_env_viz.py

```bash
python3 scripts/03_rna_spatial_simulation_with_env_viz.py
```

### Run 04_rna_spatial_simulation_env_smooth.py

```bash
python3 scripts/04_rna_spatial_simulation_env_smooth.py
```

### Run 05_rna_spatial_simulation_interactive.py

```bash
python3 scripts/05_rna_spatial_simulation_interactive.py
```

### Run 06_rna_spatial_simulation_with_statistics.py

```bash
python3 scripts/06_rna_spatial_simulation_with_statistics.py
```

### Run 07_rna_spatial_simulation_dashboard_bokeh.py

```bash
python3 -m bokeh serve --show scripts/07_rna_spatial_simulation_dashboard_bokeh.py
```

For scripts with bokeh suffix in the name:
```bash
python3 -m bokeh serve --show scripts/10_rna_spatial_simulation_dashboard_bokeh.py
```

## Screenshot

![RNA Grid Simulation Heatmap](images/heatmap_example.png)

## License

MIT License
