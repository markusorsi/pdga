"""
Configuration settings for the Peptide Design Genetic Algorithm (PDGA).
Update the values below to change the behavior of the algorithm.
"""

CONFIG = {
    # Query settings
    'query': 'NC(CCSC)C(=O)N(CC1=CC=C(Cl)C=C1(Cl))CC(=O)N(CCCCCCCC[NH3+])CC(=O)NC(CC[NH3+])C(=O)N(CC1=CC=CC=N1)CC(=O)N(CCCCO)CC(=O)O',
    'query_format': 'smiles',

    # Population settings
    'pop_size': 50,
    'pop_selection': 10,

    # Genetic operators parameters
    'mutation_ratio': 0.8,
    'cutoff': 0.5,
    'fitness_function': 'atompair',
    'selection_strategy': 'greedy',
    'crossover_method': 'skip',

    # Algorithm runtime settings
    'n_iterations': 100,
    'seed': 42,

    # Output settings
    'run_id': 'ccr8',
    'maximize': False,
}