from abc import ABC, abstractmethod
from typing import List, Dict, Type

class SelectionOperator(ABC):
    """
    Abstract base class for selection operators in the PDGA framework.

    This class defines the interface for selection operators responsible for choosing parent individuals
    from a population based on their fitness scores. Subclasses must implement the `select` method.

    Class Attributes:
        registry (Dict[str, Type['SelectionOperator']]): A dictionary mapping selection operator names
            (in lowercase) to their corresponding SelectionOperator subclass.
        name (str): The unique name identifier for the selection operator. Must be defined in subclasses.
    """

    registry: Dict[str, Type['SelectionOperator']] = {}
    name: str = None

    def __init_subclass__(cls, **kwargs):
        """
        Automatically register subclasses of SelectionOperator.

        This method ensures that any subclass of SelectionOperator defines a non-None `name` attribute.
        Upon subclass initialization, the subclass is automatically registered in the `registry` using
        its lowercase name as the key.

        Raises:
            ValueError: If a subclass does not define a `name` attribute.
        """
        super().__init_subclass__(**kwargs)
        if cls.name is None:
            raise ValueError("Subclasses of SelectionOperator must define a 'name' attribute.")
        # Automatically register the subclass using its lowercase name.
        cls.registry[cls.name.lower()] = cls

    @abstractmethod
    def select(self, fitness_scores: List[float], population: List[str], num_parents: int) -> List[str]:
        """
        Select individuals from the population based on their fitness scores.

        This abstract method must be implemented by all subclasses. It defines the logic for choosing
        a subset of the population to serve as parents for the next generation, based on the provided
        fitness scores.

        Args:
            fitness_scores (List[float]): A list of fitness scores corresponding to each individual in the population.
            population (List[str]): The current population of individuals.
            num_parents (int): The number of parent individuals to select.

        Returns:
            List[str]: A list of selected parent individuals.
        """
        pass
