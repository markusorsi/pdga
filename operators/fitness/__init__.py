from .fitness import FitnessOperator
from .fit_map4c import MAP4C
from .fit_atompair import Atompair
from .fit_mxfp import MXFP

def get_fitness_function(name: str):
    """
    Retrieve the fitness function based on the provided name.
    
    Parameters:
        name (str): The name of the fitness function (e.g., 'map4c', 'atompair', 'mxfp').
        
    Returns:
        tuple: (process_query, process, fitness)
    """
    try:
        fitness_class = FitnessOperator.registry[name.lower()]
    except KeyError:
        raise ValueError(f"Fitness function '{name}' is not defined.")
    return fitness_class()

def calculate_fitness(population, translation_dict, process_individual, fitness_function, query_processed):
    """
    Calculate fitness scores for the given population.
    
    This function replaces the evaluate_fitness method in PDGA. It processes each individual,
    computes its fitness relative to the processed query, and returns a list of fitness scores.
    
    Parameters:
        population (List[str]): List of individual sequences.
        translation_dict (dict): Dictionary for translating sequences into SMILES.
        process_individual (Callable): Function to process an individual (typically provided by the fitness function).
        fitness_function (Callable): Function that calculates fitness given a processed individual and the processed query.
        query_processed: The processed query (e.g., fingerprint) that individuals are compared against.
    
    Returns:
        List[float]: Fitness scores for each individual.
    """
    fitness_scores = []
    for individual in population:
        processed = process_individual(individual, translation_dict)
        score = fitness_function(processed, query_processed)
        fitness_scores.append(score)
    return fitness_scores