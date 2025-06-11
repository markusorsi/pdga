from .selection import SelectionOperator
from .select_greedy import Greedy
from .select_elitism import Elitism
from .select_randomize import Randomize

def get_selection_method(name: str, maximize: bool = False):
    """
    Retrieve the selection method based on the provided strategy name.

    This function looks up the appropriate selection operator from the SelectionOperator registry
    using the provided name (case-insensitive), instantiates the operator, and returns its 'select' method.

    Args:
        name (str): The name of the selection strategy (e.g., 'maximize', 'minimize', or 'randomize').

    Returns:
        function: The 'select' method of the corresponding selection operator instance.

    Raises:
        ValueError: If the specified selection method is not defined in the registry.
    """
    try:
        op_class = SelectionOperator.registry[name.lower()]
    except KeyError:
        raise ValueError(f"Selection method '{name}' is not defined.")
    return op_class().select

def select(fitness_scores, population, num_parents, maximize, selection_strategy: str = "greedy"):
    """
    Apply the specified selection strategy to choose parent individuals from the population.

    This function retrieves the appropriate selection method using `get_selection_method` and applies it
    to the provided fitness scores and population to select the specified number of parents.

    Args:
        fitness_scores (List[float]): List of fitness scores corresponding to each individual in the population.
        population (List[str]): The current population of individuals.
        num_parents (int): The number of parent individuals to select.
        selection_strategy (str, optional): The selection strategy name. Defaults to "elitism".

    Returns:
        List[str]: A list of selected parent individuals.
    """
    selection_fn = get_selection_method(selection_strategy, maximize=maximize)
    return selection_fn(fitness_scores, population, num_parents)
