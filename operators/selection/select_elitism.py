from .selection import SelectionOperator
from typing import List

class Elitism(SelectionOperator):
    """
    A selection operator that applies elitism by always retaining a fixed number of top individuals.

    This operator is suitable for both minimization and maximization tasks and requires specifying
    whether to sort in ascending or descending order based on the goal.

    Attributes:
        name (str): The unique identifier for this selection operator.
        elite_count (int): Number of top individuals to always retain.
        minimize (bool): Whether the objective is minimization (True) or maximization (False).
    """
    name = 'elitism'

    def __init__(self, elite_count: int = 5, maximize: bool = True):
        """
        Initialize the Elitism operator with sorting direction and elite count.

        Args:
            elite_count (int): Number of top individuals to retain.
            minimize (bool): Whether to minimize (True) or maximize (False) the fitness score.
        """
        self.elite_count = elite_count
        self.maximize = maximize

    def select(self, fitness_scores: List[float], population: List[str], num_parents: int) -> List[str]:
        """
        Select parent individuals using elitism and fitness-based ranking.

        The top `elite_count` individuals are always selected. The remaining parents are selected
        from the rest of the population based on fitness score.

        Args:
            fitness_scores (List[float]): A list of fitness scores corresponding to each individual in the population.
            population (List[str]): The current population of individuals.
            num_parents (int): The total number of parents to select (including elites).

        Returns:
            List[str]: A list of selected parent individuals.
        """
        ranked_indices = sorted(range(len(fitness_scores)), key=lambda idx: fitness_scores[idx], reverse=self.maximize)
        elites = [population[idx] for idx in ranked_indices[:self.elite_count]]
        remaining = [population[idx] for idx in ranked_indices[self.elite_count:num_parents]]
        return elites + remaining