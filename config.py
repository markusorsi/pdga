"""
Configuration settings for the Peptide Design Genetic Algorithm (PDGA).
Update the values below to change the behavior of the algorithm.
"""

CONFIG = {
    # Query settings
    'query': 'CC(C)CCCCC(=O)N[C@@H](CC[NH3+])C(=O)N[C@@H]([C@H](C)O)C(=O)N[C@@H](CC[NH3+])C(=O)N[C@@H](CCN8)C(=O)N[C@@H](CC[NH3+])C(=O)N[C@H](CC1=CC=CC=C1)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC[NH3+])C(=O)N[C@@H](CC[NH3+])C(=O)N[C@@H]([C@H](C)O)C(=O)8',
    'query_format': 'smiles',

    # Population settings
    'pop_size': 50,
    'pop_selection': 10,

    # Genetic operators parameters
    'mutation_ratio': 0.8,
    'cutoff': 0.9,
    'fitness_function': 'atompair',
    'selection_strategy': 'greedy',
    'crossover_method': 'single_point',

    # Algorithm runtime settings
    'n_iterations': 1000,
    'seed': 42,

    # Output settings
    'run_id': 'ccr8',
    'maximize': False,
}