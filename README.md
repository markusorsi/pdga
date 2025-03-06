# PDGA - Peptide Design Genetic Algorithm

PDGA is a modular genetic algorithm framework designed to generate peptide sequences with desired properties. It integrates modular operators for fitness evaluation, selection, crossover, and mutation. In addition, PDGA uses a flexible configuration system, a building block manager for processing chemical building blocks, and a dedicated results handler to collect and export experiment hits.

## Features

- **Modular Design:** Easily swap out fitness, selection, crossover, and mutation operators.
- **Building Block Management:** Use CSV files to define building blocks, N-caps, branches. A dedicated `BuildingBlockManager` handles data loading and sequence translation.
- **Query Handling:** `QueryHandler` processes the input query based on its format (e.g., SMILES or building block sequence) into a comparable format.
- **Results Handling:** `ResultsHandler` collects hit sequences (those meeting a fitness cutoff) and exports them to a CSV file with configurable sorting.
- **Configurable via CLI:** Experiment parameters can be provided via a configuration file as well as command-line arguments, enabling easy switching between experiments.

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/yourusername/pdga.git
   cd pdga
   ```

2. **Create a Virtual Environment and Install Dependencies:**

   ```bash
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```

3. **Prepare Your Data:**

   Place your CSV files containing building block data in the `bblocks/` directory:
   - `bb_monomers.csv`
   - `bb_ncaps.csv`
   - `branches.csv`
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
- `--seed`: Random seed for reproducibility.
- `--sort_ascending`: Sort results in ascending order (if provided; otherwise, descending).

## File Structure

```
pdga/
├── bblocks/
│   ├── bb_monomers.csv
│   ├── bb_ncaps.csv
│   ├── bb_additional.csv
│   ├── branches.csv
│   ├── bbmanager.py         # BuildingBlockManager implementation.
│   └── __init__.py
├── config.py                # Default configuration file.
├── main.py                  # Entry point script.
├── operators/
│   ├── crossover/           # Crossover operators and helper functions.
│   ├── fitness/             # Fitness functions and abstract base class.
│   ├── mutation/            # Mutation operators.
│   └── selection/           # Selection operators.
├── pdga.py                  # PDGA class implementation.
├── requirements.txt         # Python package dependencies.
├── README.md                # This file.
└── utils.py                 # Utility functions (e.g., ResultsHandler, QueryHandler).
```

## Configuration File Example

Below is a sample `config.py` file:

```python
CONFIG = {
    # Query settings
    'query': 'CC(C)CCCCC(=O)N[C@@H](CCN)C(=O)...',  # truncated SMILES for target molecule
    'query_format': 'smiles',

    # Population settings
    'pop_size': 50,
    'pop_selection': 10,

    # Genetic operators parameters
    'mutation_ratio': 0.5,
    'cutoff': 0.5,
    'fitness_function': 'atompair',  # Options: 'map4c', 'atompair', 'mxfp'
    'selection_strategy': 'maximize',
    'crossover_method': 'single_point',

    # Algorithm runtime settings
    'n_iterations': 300,
    'seed': 42,

    # Output settings
    'run_id': 'experiment_01',
    'sort_ascending': False  # False sorts results in descending order (higher fitness first)
}
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or feedback, please contact [markusorsi@icloud.com](mailto:markusorsi@icloud.com).