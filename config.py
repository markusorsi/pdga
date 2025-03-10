"""
Configuration settings for the Peptide Design Genetic Algorithm (PDGA).
Update the values below to change the behavior of the algorithm.
"""

CONFIG = {
    # Query settings
    'query': 'SEQVENCE',
    'query_format': 'sequence',

    # Population settings
    'pop_size': 50,
    'pop_selection': 10,

    # Genetic operators parameters
    'mutation_ratio': 0.5,
    'cutoff': 300,
    'fitness_function': 'mxfp',
    'selection_strategy': 'minimize',
    'crossover_method': 'single_point',

    # Algorithm runtime settings
    'n_iterations': 1000,
    'seed': 42,

    # Output settings
    'run_id': 'test_run',
    'maximize': False,
}