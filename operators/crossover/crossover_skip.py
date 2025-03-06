from typing import Tuple
from .crossover import CrossoverOperator

class Skip(CrossoverOperator):
    name = 'skip'

    def crossover(self, parent1: str, parent2: str) -> Tuple[str, str]:
        """
        Skip crossover by returning the two parent sequences without any modification.

        Parameters:
            parent1 (str): The first parent sequence.
            parent2 (str): The second parent sequence.

        Returns:
            Tuple[str, str]: The unchanged parent sequences.
        """
        return parent1, parent2
