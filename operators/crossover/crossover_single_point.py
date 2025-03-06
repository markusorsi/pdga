import random
from typing import Tuple
from .crossover import CrossoverOperator

class SinglePoint(CrossoverOperator):
    name = 'single_point'

    def crossover(self, parent1: str, parent2: str) -> Tuple[str, str]:
        """
        Perform single-point crossover between two parent sequences.
        
        Parameters:
            parent1 (str): The first parent sequence.
            parent2 (str): The second parent sequence.
        
        Returns:
            Tuple[str, str]: Two offspring produced from the parents.
        """
        parent1_parts = parent1.split('-')
        parent2_parts = parent2.split('-')

        # If either parent has fewer than 2 parts, crossover is not possible.
        if len(parent1_parts) < 2 or len(parent2_parts) < 2:
            return parent1, parent2

        # Choose a random crossover point.
        crossover_point = random.randint(1, len(parent1_parts) - 1)
        offspring1 = '-'.join(parent1_parts[:crossover_point] + parent2_parts[crossover_point:])
        offspring2 = '-'.join(parent2_parts[:crossover_point] + parent1_parts[crossover_point:])

        # Sanitize the offspring using the base class's method.
        offspring1 = self.sanitize_offspring(offspring1)
        offspring2 = self.sanitize_offspring(offspring2)

        return offspring1, offspring2
