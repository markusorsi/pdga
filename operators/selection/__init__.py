from .selection import SelectionOperator
from .select_maximize import Maximize
from .select_minimize import Minimize
from .select_randomize import Randomize

def get_selection_method(name: str):
    """
    Retrieve the selection method based on the provided name.

    Parameters:
        name (str): The name of the selection strategy (e.g., 'maximize' or 'minimize').

    Returns:
        function: The 'select' method of the corresponding selection operator instance.
    """
    try:
        op_class = SelectionOperator.registry[name.lower()]
    except KeyError:
        raise ValueError(f"Selection method '{name}' is not defined.")
    return op_class().select

def select(fitness_scores, population, num_parents, selection_strategy: str = "maximize"):
    """
    Apply selection to the given population using the specified strategy.

    Parameters:
        fitness_scores (List[float]): List of fitness scores.
        population (List[str]): The current population of individuals.
        num_parents (int): The number of parents to select.
        selection_strategy (str): Strategy name (e.g., 'maximize' or 'minimize').

    Returns:
        List[str]: A list of selected parent individuals.
    """
    selection_fn = get_selection_method(selection_strategy)
    return selection_fn(fitness_scores, population, num_parents)