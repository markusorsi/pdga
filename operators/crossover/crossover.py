from abc import ABC, abstractmethod
from random import randint
import re

class CrossoverOperator(ABC):
    registry = {}
    name: str = None  # Subclasses must define a unique name.

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name is None:
            raise ValueError('Subclasses of CrossoverOperator must define a "name" attribute.')
        # Register the subclass using its lowercase name.
        cls.registry[cls.name.lower()] = cls

    @abstractmethod
    def crossover(self, parent1: str, parent2: str) -> tuple:
        """
        Perform crossover between two parent sequences and return two offspring.
        
        Parameters:
            parent1 (str): The first parent sequence.
            parent2 (str): The second parent sequence.
        
        Returns:
            tuple: A tuple containing two offspring sequences.
        """
        pass

    def sanitize_offspring(self, offspring: str) -> str:
        """
        Sanitize an offspring sequence to enforce feasibility rules.

        This method applies several modifications to the input sequence to ensure that it meets chemical feasibility constraints:

        1. Disulfide Bonds ('s'):
            - If the number of 's' tokens is odd, remove one occurrence. The removal targets an 's'
            token that appears either at the beginning (as 's-'), at the end (as '-s'), or in the middle (as '-s').
        2. Branch Tokens ('b'):
            - If more than one branch token is present remove all but one occurrence.
        3. Cyclization ('c'):
            - If a cyclization token ('c') is present and there is at least one branch token ('b') or an N-cap token ('T'),
            remove the cyclization token.

        Parameters:
            offspring (str): The offspring sequence to sanitize.

        Returns:
            str: The sanitized offspring sequence with the above rules applied.
        """
        num_c = len(re.findall(r'c', offspring))
        num_b = len(re.findall(r'b', offspring))
        num_s = len(re.findall(r's', offspring))
        num_T = len(re.findall(r'T', offspring))

        # Handles cases for a single 's' token by removing one occurrence
        # from the beginning, end, or middle of the sequence.
        if num_s % 2 != 0:
            offspring = re.sub(r'(^s-)|(-s$)|(-s)', '', offspring, count=1)
        
        # If more than one branch token is present, remove the extra occurrences.
        if num_b > 1:
            offspring = re.sub(r'\bb\d{3}\b-', '', offspring, count=num_b-1)
        
        # If a cyclization token is present and there is at least one branch or N-cap token,
        # remove the cyclization token.
        if num_c > 0 and (num_b > 0 or num_T > 0):
            offspring = re.sub(r'c-', '', offspring)
        
        return offspring
