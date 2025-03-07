<img src="https://img.shields.io/badge/Python-3.10-blue?style=flat-square"/> <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square"/>

# Peptide Design Genetic Algorithm (PDGA)

PDGA is a modular genetic algorithm framework designed to generate peptide sequences with desired building blocks and properties. It integrates modular operators for fitness evaluation, selection, crossover, and mutation. In addition, PDGA uses a flexible configuration system, a building block manager for processing chemical building blocks, and a dedicated results handler to collect and export experimental hits.

## Features

- **Modular Design:** Easily swap out building blocks, fitness, selection, crossover, and mutation operators.
- **Simple Implementation:** `PDGA` is a simple implementation that can be easily adapted to different use cases for peptide design.
- **Building Block Management:** `BuildingBlockManager` handles data loading and sequence translation of building blocks 
- **Results Handling:** `ResultsHandler` collects hit sequences (those meeting a fitness cutoff) and exports them to a CSV file with configurable sorting.
- **Configurable via CLI:** Experiment parameters can be provided via a configuration file as well as command-line arguments, enabling easy switching between experiments.

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/markusorsi/pdga.git
   cd pdga
   ```

2. **Create a Virtual Environment and Install Dependencies:**

   ```bash
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```

3. **Prepare Your Data:**

   Replace of modify the CSV files containing building blocks in the `bblocks/` directory:
   - `bb_monomers.csv`
   - `bb_ncaps.csv`
   - `bb_branches.csv`
   - `bb_additional.csv`

## Usage

### Running PDGA

You can run PDGA from the command line. The main script accepts various command-line arguments to override configuration defaults.

For example:

```bash
./main.py --config config.py --pop_size 100 --n_iterations 500 --run_id experiment_01
```

The `--config` parameter allows you to specify a configuration file containing default settings. All parameters (such as query, population size, mutation ratio, fitness function, selection strategy, crossover method, etc.) can be overridden via the terminal.

### Command-Line Arguments

- `--config`: Path to the configuration file (default: `config.py`).
- `--query`: Target query used for fitness evaluation.
- `--query_format`: Format of the query (e.g., `smiles`).
- `--pop_size`: Population size.
- `--pop_selection`: Number of individuals selected as parents.
- `--mutation_ratio`: Ratio of individuals to be mutated.
- `--cutoff`: Fitness cutoff to record hits.
- `--fitness_function`: Identifier for the fitness function to use (e.g., `map4c`, `atompair`, `mxfp`).
- `--selection_strategy`: Identifier for the selection strategy (e.g., `maximize`).
- `--crossover_method`: Identifier for the crossover method (e.g., `single_point`, `skip`).
- `--n_iterations`: Number of generations to run.
- `--run_id`: Identifier for the run (used in output filenames).
- `--maximize`: The goal of the optimization study.
- `--seed`: Random seed for reproducibility.

## File Structure
```
├── bblocks
│   ├── bb_additional.csv     # Additional building blocks.
│   ├── bb_branches.csv       # List of branch building blocks.
│   ├── bbmanager.py          # BuildingBlockManager implementation.
│   ├── bb_monomers.csv       # List of monomers concatenable to form peptides.
│   └── bb_ncaps.csv          # List of N-terminal caps.
├── config.py                 # Default configuration file.
├── logs                      # Log directory for experiments.
├── main.py                   # Entry point script.
├── operators
│   ├── crossover             # Crossover operators.
│   ├── fitness               # Fitness functions.
│   ├── mutation              # Mutation operators.
│   └── selection             # Selection operators.
├── pdga.py                   # Main PDGA implementation.
├── README.md                 # This file.  
├── requirements.txt          # Python dependencies.
├── results                   # Output directory for experiments.
└── utils.py                  # ResultsHandler implementation.
```

## Configuration File Example

Below is a sample `config.py` file:

```python
CONFIG = {
    # Query settings
    'query': 'SEQVENCE',
    'query_format': 'sequence',

    # Population settings
    'pop_size': 50,
    'pop_selection': 10,

    # Genetic operators parameters
    'mutation_ratio': 0.5,
    'cutoff': 0.5,
    'fitness_function': 'atompair',
    'selection_strategy': 'maximize',
    'crossover_method': 'single_point',

    # Algorithm runtime settings
    'n_iterations': 1000,
    'seed': 42,

    # Output settings
    'run_id': 'test_run',
    'maximize': True,
}
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or feedback, please contact [markusorsi@icloud.com](mailto:markusorsi@icloud.com).

## Citation

```
@article{orsi2025navigating,
  title = {Navigating a 1E+60 Chemical Space of Peptide/Peptoid Oligomers},
  author = {Orsi, Markus and Reymond, Jean-Louis},
  journal = {Molecular Informatics},
  volume = {44},
  number = {1},
  pages = {e202400186},
  year = {2025},
  publisher = {Wiley Online Library}
}
```