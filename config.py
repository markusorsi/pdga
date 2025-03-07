"""
Configuration settings for the Peptide Design Genetic Algorithm (PDGA).
Update the values below to change the behavior of the algorithm.
"""

CONFIG = {
    # Query settings
    'query': 'CC(C)CCCCC(=O)N[C@@H](CCN)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CCN)C(=O)N[C@H]1CCNC(=O)[C@@H](NC(=O)[C@H](CCN)NC(=O)[C@H](CCN)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](Cc2ccccc2)NC(=O)[C@H](CCN)NC1=O)[C@@H](C)O',
    'query_format': 'smiles',

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
    'n_iterations': 10000,
    'seed': 42,

    # Output settings
    'run_id': 'polymyxin',
    'sort_ascending': False,
}
