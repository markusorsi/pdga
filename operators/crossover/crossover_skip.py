from typing import Tuple
from .crossover import CrossoverOperator

class Skip(CrossoverOperator):
    """
    A crossover operator that performs no crossover.

    This operator simply returns the parent sequences unmodified, effectively skipping
    the crossover process.
    
    Attributes:
        name (str): Unique identifier for this crossover operator.
    """
    name = 'skip'

    def crossover(self, parent1: str, parent2: str) -> Tuple[str, str]:
        """
        Skip crossover by returning the two parent sequences without any modification.

        Args:
            parent1 (str): The first parent sequence.
            parent2 (str): The second parent sequence.

        Returns:
            Tuple[str, str]: A tuple containing the unchanged parent sequences.
        """
        return parent1, parent2
