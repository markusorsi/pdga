from .selection import SelectionOperator
from typing import List
import random

class Randomize(SelectionOperator):
    """
    A selection operator that randomly selects parent individuals from the population.

    This operator disregards fitness scores and chooses parents at random.
    
    Attributes:
        name (str): The unique identifier for this selection operator.
    """
    name = 'randomize'

    def select(self, fitness_scores: List[float], population: List[str], num_parents: int) -> List[str]:
        """
        Randomly select parent individuals from the population.

        This method randomly samples a specified number of individuals from the population,
        ignoring the provided fitness scores.

        Args:
            fitness_scores (List[float]): A list of fitness scores (unused in this implementation).
            population (List[str]): The current population of individuals.
            num_parents (int): The number of parent individuals to select.

        Returns:
            List[str]: A list of randomly selected parent individuals.
        """
        selected_parents = random.sample(population, num_parents)
        return selected_parents
