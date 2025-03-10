from abc import ABC, abstractmethod
from typing import Any, Dict, Type

class FitnessOperator(ABC):
    """
    Abstract base class for fitness functions.

    Subclasses must implement the following methods:
      - process_query: Convert the query into a processed representation (e.g., a fingerprint).
      - process: Convert an individual sequence into a comparable format using a translation dictionary.
      - fitness: Calculate a fitness score for an individual given the processed query.

    Class Attributes:
        registry (Dict[str, Type['FitnessOperator']]): A dictionary mapping fitness operator names 
            (in lowercase) to their corresponding FitnessOperator subclasses.
        name (str): The unique name identifier for the fitness operator. Must be defined in subclasses.
    """
    registry: Dict[str, Type['FitnessOperator']] = {}
    name: str = None

    def __init_subclass__(cls, **kwargs):
        """
        Automatically register subclasses of FitnessOperator.

        Ensures that any subclass of FitnessOperator defines a non-None `name` attribute.
        Upon subclass initialization, the subclass is automatically registered in the `registry`
        using its lowercase name as the key.

        Raises:
            ValueError: If a subclass does not define a `name` attribute.
        """
        super().__init_subclass__(**kwargs)
        if cls.name is None:
            raise ValueError("Subclasses of FitnessOperator must define a 'name' attribute.")
        # Automatically register the subclass using its name in lowercase.
        cls.registry[cls.name.lower()] = cls

    @abstractmethod
    def process_query(self, query: str, query_format: str) -> Any:
        """
        Process the query and return a processed representation.

        This method should transform the raw query into a format suitable for fitness evaluation,
        such as converting a SMILES string into a molecular fingerprint.

        Args:
            query (str): The raw query input.
            query_format (str): The format of the query (e.g., 'smiles').

        Returns:
            Any: A processed representation of the query.
        """
        pass

    @abstractmethod
    def process(self, sequence: str, translation_dict: Dict) -> Any:
        """
        Process a sequence using a translation dictionary.

        Converts the sequence into a format that can be compared with the processed query.
        This may involve translating tokens in the sequence using the provided translation dictionary.

        Args:
            sequence (str): The sequence to process.
            translation_dict (Dict): A dictionary mapping sequence tokens to their corresponding representations.

        Returns:
            Any: A processed representation of the sequence.
        """
        pass

    @abstractmethod
    def fitness(self, individual: Any, query: Any) -> float:
        """
        Calculate and return the fitness score for an individual.

        The fitness score is computed based on the processed individual and the processed query.
        Higher (or lower, depending on the goal) fitness scores indicate a better match.

        Args:
            individual (Any): The processed representation of a sequence.
            query (Any): The processed representation of the query.

        Returns:
            float: The fitness score of the individual.
        """
        pass
