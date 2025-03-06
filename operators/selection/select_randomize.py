from .selection import SelectionOperator
from typing import List
import random

class Randomize(SelectionOperator):
    name = 'randomize'

    def select(self, fitness_scores: List[float], population: List[str], num_parents: int) -> List[str]:
        """
        Select individuals by choosing randomly.
        """
        selected_parents = random.sample(population, num_parents)
        return selected_parents