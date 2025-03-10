from .fitness import FitnessOperator
from .fit_map4c import MAP4C
from .fit_atompair import Atompair
from .fit_mxfp import MXFP

def get_fitness_function(name: str):
    """
    Retrieve and instantiate the fitness function based on the provided name.

    This function looks up the appropriate FitnessOperator subclass in the registry using the
    provided name (case-insensitive) and returns an instance of the corresponding fitness operator.

    Args:
        name (str): The name of the fitness function (e.g., 'map4c', 'atompair', 'mxfp').

    Returns:
        FitnessOperator: An instance of the corresponding fitness operator.

    Raises:
        ValueError: If the specified fitness function is not defined in the registry.
    """
    try:
        fitness_class = FitnessOperator.registry[name.lower()]
    except KeyError:
        raise ValueError(f"Fitness function '{name}' is not defined.")
    return fitness_class()

def calculate_fitness(population, translation_dict, process_individual, fitness_function, query_processed):
    """
    Calculate fitness scores for a population of individuals.

    This function processes each individual sequence using the provided processing function and
    translation dictionary, then computes the fitness score relative to the processed query using
    the given fitness function. It returns a list of fitness scores corresponding to each individual.

    Args:
        population (List[str]): List of individual sequences.
        translation_dict (dict): Dictionary for translating sequences into SMILES representations.
        process_individual (Callable): Function to process an individual sequence (typically provided by the fitness operator).
        fitness_function (Callable): Function that calculates fitness given a processed individual and the processed query.
        query_processed: The processed query (e.g., a molecular fingerprint) against which individuals are compared.

    Returns:
        List[float]: Fitness scores for each individual in the population.
    """
    fitness_scores = []
    for individual in population:
        processed = process_individual(individual, translation_dict)
        score = fitness_function(processed, query_processed)
        fitness_scores.append(score)
    return fitness_scores
