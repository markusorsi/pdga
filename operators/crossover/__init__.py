from .crossover import CrossoverOperator
from .crossover_single_point import SinglePoint
from .crossover_skip import Skip

def get_crossover_method(name: str):
    """
    Retrieve the crossover method based on the provided strategy name.

    This function looks up the appropriate CrossoverOperator subclass in the registry using the
    provided name (case-insensitive), instantiates the operator, and returns its crossover function.

    Args:
        name (str): The name of the crossover strategy (e.g., 'single_point' or 'skip').

    Returns:
        function: The crossover method of the corresponding crossover operator instance.

    Raises:
        ValueError: If the specified crossover method is not defined.
    """
    try:
        op_class = CrossoverOperator.registry[name.lower()]
    except KeyError:
        raise ValueError(f"Crossover method '{name}' is not defined.")
    return op_class().crossover
