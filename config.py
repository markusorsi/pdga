"""
Configuration settings for the Peptide Design Genetic Algorithm (PDGA).
Update the values below to change the behavior of the algorithm.
"""

CONFIG = {
    # Query settings
    'query': 'KLQwWkLKPlKL',
    'query_format': 'sequence',

    # Population settings
    'pop_size': 50,
    'pop_selection': 10,

    # Genetic operators parameters
    'mutation_ratio': 0.5,
    'cutoff': 250,
    'fitness_function': 'atompair',
    'selection_strategy': 'maximize',
    'crossover_method': 'single_point',

    # Algorithm runtime settings
    'n_iterations': 100,
    'seed': 42,

    # Output settings
    'run_id': 'polymyxin',
    'maximize': True,
}