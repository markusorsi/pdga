from abc import ABC, abstractmethod
from random import randint
import re

class CrossoverOperator(ABC):
    """
    Abstract base class for crossover operators in the PDGA framework.

    This class defines the interface for performing crossover between two parent sequences
    to produce offspring. Subclasses must implement the `crossover` method and define a unique
    `name` attribute. Each subclass is automatically registered in the `registry` dictionary
    using its lowercase name.

    Class Attributes:
        registry (dict): A dictionary mapping crossover operator names (in lowercase) to their corresponding subclasses.
        name (str): Unique identifier for the crossover operator (must be defined in subclasses).
    """
    registry = {}
    name: str = None  # Subclasses must define a unique name.

    def __init_subclass__(cls, **kwargs):
        """
        Automatically register subclasses of CrossoverOperator.

        This method checks that the subclass defines a non-None `name` attribute and registers it
        in the `registry` dictionary using its lowercase name as the key.

        Raises:
            ValueError: If a subclass does not define a `name` attribute.
        """
        super().__init_subclass__(**kwargs)
        if cls.name is None:
            raise ValueError('Subclasses of CrossoverOperator must define a "name" attribute.')
        # Register the subclass using its lowercase name.
        cls.registry[cls.name.lower()] = cls

    @abstractmethod
    def crossover(self, parent1: str, parent2: str) -> tuple:
        """
        Perform crossover between two parent sequences and return two offspring.

        This abstract method must be implemented by subclasses to define the specific crossover strategy.

        Args:
            parent1 (str): The first parent sequence.
            parent2 (str): The second parent sequence.

        Returns:
            tuple: A tuple containing two offspring sequences.
        """
        pass

    def sanitize_offspring(self, offspring: str) -> str:
        """
        Sanitize an offspring sequence to enforce chemical feasibility constraints.

        This method applies several modifications to the input sequence to ensure that it meets feasibility rules:
        
        1. Disulfide Bonds ('s'):
           - If the number of 's' tokens is odd, remove one occurrence. The removal targets an 's' token that appears
             either at the beginning (as 's-'), at the end (as '-s'), or in the middle (as '-s').
        2. Branch Tokens ('b'):
           - If more than one branch token is present, remove the extra occurrences, keeping only one.
        3. Cyclization ('c'):
           - If a cyclization token ('c') is present and there is at least one branch token ('b') or an N-cap token ('T'),
             remove the cyclization token.

        Args:
            offspring (str): The offspring sequence to sanitize.

        Returns:
            str: The sanitized offspring sequence with the feasibility rules applied.
        """
        num_c = len(re.findall(r'c', offspring))
        num_b = len(re.findall(r'b', offspring))
        num_s = len(re.findall(r's', offspring))
        num_T = len(re.findall(r'T', offspring))

        # Handles cases for an odd number of 's' tokens by removing one occurrence from
        # the beginning, end, or middle of the sequence.
        if num_s % 2 != 0:
            offspring = re.sub(r'(^s-)|(-s$)|(-s)', '', offspring, count=1)
        
        # If more than one branch token is present, remove the extra occurrences.
        if num_b > 1:
            offspring = re.sub(r'\bb\d{3}\b-', '', offspring, count=num_b-1)
        
        # If a cyclization token is present and there is at least one branch or N-cap token, remove the cyclization token.
        if num_c > 0 and (num_b > 0 or num_T > 0):
            offspring = re.sub(r'c-', '', offspring)
        
        return offspring
