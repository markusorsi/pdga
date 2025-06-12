<img src="https://img.shields.io/badge/Python-3.12.7-blue?style=flat-square"/> <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square"/>

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
   python -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```

   _Tip: The python version used in this project is displayed in the badge above. Make sure to use the same version to avoid compatibility issues._

3. **Edit Your Building Blocks:**

PDGA constructs peptides from modular building blocks defined in CSV files located in the `bblocks/` directory. These files include:

- `bb_monomers.csv`
- `bb_ncaps.csv`
- `bb_branches.csv`  
- `bb_additional.csv`

Each file can be modified or extended to customize the design space.

#### `bb_monomers.csv`

This file contains the main building blocks (standard amino acids and analogs). Each entry includes:

- A unique ID in the format `BBXXX` (e.g., `BB901`, `BB902`, etc.)
- A valid SMILES string
- An **optional** human-readable label

Example:
```
ID       SMILES                             Description
BB901    NC(C(C)C)C(=O)                     Valine
BB902    NC(CC(C)C)C(=O)                    Leucine
BB903    NC(C(C)CC)C(=O)                    Isoleucine
```

- SMILES must **start with the N-terminus** and **end with the C-terminus** (without the final `O`).
- The algorithm automatically adds the terminal `O` during final sequence generation.
- The ID format (`BBXXX`) must be preserved; the number after `BB` can be customized.

---

#### `bb_ncaps.csv`

This file defines N-terminal caps. Each entry:

- Should use the ID format `TXXX` (e.g., `T001`, `T002`)
- Must **end in a C-terminus** (e.g., a carboxyl group)

---

#### `bb_branches.csv`

This file defines branching building blocks used for side-chain cyclization. Each entry:

- Uses an ID in the format `bXXX` (e.g., `b001`)
- Must contain a SMILES string with a **side-chain cyclization point**
- Must mark the **cyclizing atom with an `8`** to indicate where the algorithm should connect

---

## Usage

### Running PDGA

You can run PDGA from the command line. The main script accepts various command-line arguments to override configuration defaults.

For example:

```bash
python main.py
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
- `--selection_strategy`: Identifier for the selection strategy (e.g., `greedy`).
- `--crossover_method`: Identifier for the crossover method (e.g., `single_point`, `skip`).
- `--n_iterations`: Number of generations to run.
- `--run_id`: Identifier for the run (used in output filenames).
- `--maximize`: The goal of the optimization study.
- `--seed`: Random seed for reproducibility.

---

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

---

## Configuration File Example

Below is a sample `config.py` file:

```python
CONFIG = {
    # === Query Settings ===
    
    # Target query used for fitness evaluation. Should match the selected fitness function's format.
    'query': 'CC(C)CCCCC(=O)N[C@@H](CC[NH3+])C(=O)N[C@@H]([C@H](C)O)C(=O)N[C@@H](CC[NH3+])C(=O)N[C@@H](CCN8)C(=O)N[C@@H](CC[NH3+])C(=O)N[C@H](CC1=CC=CC=C1)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC[NH3+])C(=O)N[C@@H](CC[NH3+])C(=O)N[C@@H]([C@H](C)O)C(=O)8',  
    'query_format': 'smiles',  # Options: 'smiles', 'sequence', etc.

    # === Population Settings ===

    'pop_size': 50,           # Total number of individuals maintained in the population.
    'pop_selection': 10,      # Number of individuals selected each generation for mutation/crossover.

    # === Genetic Operator Parameters ===

    'mutation_ratio': 0.8,    # Probability that an offspring undergoes mutation (0 to 1).
    'cutoff': 0.5,            # Fitness threshold to record a sequence as a "hit" in the results file.

    # Fitness functions to choose from:
    # 'map4c'     - MAP4 fingerprint + Jaccard distance (minimize)
    # 'mxfp'      - MXFP fingerprint + Manhattan distance (minimize)
    # 'atompair'  - Atom pair fingerprint + Jaccard distance (minimize)
    'fitness_function': 'atompair',

    # Selection strategies available:
    # 'greedy'    - Top-scoring individuals selected deterministically
    # 'randomize' - Random selection from the population 
    'selection_strategy': 'greedy',

    # Crossover methods to choose from:
    # 'single_point' - Simple one-point crossover
    # 'skip'         - No crossover
    'crossover_method': 'single_point',

    # === Runtime Settings ===

    'n_iterations': 1000,     # Number of generations to evolve
    'seed': 42,               # Random seed for reproducibility

    # === Output Settings ===

    'run_id': 'run',          # Custom label for the run; used in output filenames
    'maximize': False         # Set to True if higher fitness = better; False to minimize
}
```

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contact

For questions or feedback, please contact [markusorsi@icloud.com](mailto:markusorsi@icloud.com).

---

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
