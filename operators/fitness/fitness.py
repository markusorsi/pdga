from abc import ABC, abstractmethod
from typing import Any, Dict, Type

class FitnessOperator(ABC):
    """
    Abstract base class for fitness functions.
    Subclasses must implement the following methods:
      - process_query: Convert the query into a format (e.g. a fingerprint).
      - process: Convert an individual (sequence) into a comparable format.
      - fitness: Calculate a fitness score given an individual and a processed query.
    """
    registry: Dict[str, Type['FitnessOperator']] = {}
    name: str = None

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name is None:
            raise ValueError("Subclasses of FitnessOperator must define a 'name' attribute.")
        # Automatically register the subclass using its name in lowercase.
        cls.registry[cls.name.lower()] = cls

    @abstractmethod
    def process_query(self, query: str) -> Any:
        """
        Process the query and return a processed representation.
        """
        pass

    @abstractmethod
    def process(self, sequence: str, translation_dict: Dict) -> Any:
        """
        Process a peptide sequence using a translation dictionary.
        """
        pass

    @abstractmethod
    def fitness(self, individual: Any, query: Any) -> float:
        """
        Calculate and return the fitness score given a processed individual and query.
        """
        pass
