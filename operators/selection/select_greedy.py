from .selection import SelectionOperator
from typing import List

class Greedy(SelectionOperator):
    """
    A greedy selection operator that selects individuals with the highest fitness scores.

    This operator is used in maximization problems, where a higher fitness score indicates a better solution.
    
    Attributes:
        name (str): The unique identifier for this selection operator.
    """
    name = 'greedy'

    def __init__(self, maximize: bool = False):
        """
        Initialize the Greedy operator with sorting direction.

        Args:
            maximize (bool): Whether to maximize (True) or minimize (False) the fitness score.
        """
        self.maximize = maximize

    def select(self, fitness_scores: List[float], population: List[str], num_parents: int) -> List[str]:
        """
        Select parent individuals based on the highest fitness scores.

        This method ranks the population by fitness score in descending order (highest scores first)
        and selects the top `num_parents` individuals.

        Args:
            fitness_scores (List[float]): A list of fitness scores corresponding to each individual in the population.
            population (List[str]): The current population of individuals.
            num_parents (int): The number of parent individuals to select.

        Returns:
            List[str]: A list of individuals with the highest fitness scores, selected as parents.
        """
        ranked_indices = sorted(range(len(fitness_scores)), key=lambda idx: fitness_scores[idx], reverse=self.maximize)
        selected_parents = [population[idx] for idx in ranked_indices[:num_parents]]
        return selected_parents
