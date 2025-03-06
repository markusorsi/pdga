from .selection import SelectionOperator
from typing import List

class Minimize(SelectionOperator):
    name = 'minimize'

    def select(self, fitness_scores: List[float], population: List[str], num_parents: int) -> List[str]:
        """
        Select individuals by choosing the ones with the lowest fitness scores.
        """
        ranked_indices = sorted(range(len(fitness_scores)), key=lambda idx: fitness_scores[idx])
        selected_parents = [population[idx] for idx in ranked_indices[:num_parents]]
        return selected_parents