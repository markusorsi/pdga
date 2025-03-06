from abc import ABC, abstractmethod
from typing import List,Dict, Type

class SelectionOperator(ABC):

    registry: Dict[str, Type['SelectionOperator']] = {}
    name: str = None

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name is None:
            raise ValueError("Subclasses of SelectionOperator must define a 'name' attribute.")
        # Automatically register the subclass using its name in lowercase.
        cls.registry[cls.name.lower()] = cls

    @abstractmethod
    def select(self, fitness_scores: List[float], population: List[str], num_parents: int) -> List[str]:
        """
        Select individuals from the population based on fitness scores.
        
        Parameters:
            fitness_scores (List[float]): The list of fitness scores.
            population (List[str]): The current population of individuals.
            num_parents (int): The number of parents to select.
        
        Returns:
            List[str]: A list of selected parent individuals.
        """
        pass