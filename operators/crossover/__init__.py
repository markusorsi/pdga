from .crossover import CrossoverOperator
from .crossover_single_point import SinglePoint
from .crossover_skip import Skip

def get_crossover_method(name: str):
    """
    Retrieve the crossover method based on the provided name.
    
    Parameters:
        name (str): The name of the crossover strategy (e.g., 'single_point' or 'skip').
        
    Returns:
        function: The crossover function from the corresponding operator instance.
    """
    try:
        op_class = CrossoverOperator.registry[name.lower()]
    except KeyError:
        raise ValueError(f'Crossover method \'{name}\' is not defined.')
    return op_class().crossover