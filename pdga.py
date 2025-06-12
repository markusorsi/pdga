"""
This module implements the Peptide Design Genetic Algorithm (PDGA) for designing peptide sequences.

Classes:
    PDGA: Implements the genetic algorithm for peptide design.
"""

import time
import random
from typing import List

from bblocks.bbmanager import BuildingBlockManager
from operators.fitness import get_fitness_function
from operators.selection import get_selection_method
from operators.crossover import get_crossover_method
from operators.mutation import mutate
from utils import ResultsHandler

from tqdm import tqdm

class PDGA:
    """
    Peptide Design Genetic Algorithm (PDGA)

    This class implements a generative genetic algorithm to design peptide sequences.
    It uses modular operators for fitness evaluation, selection, crossover, and mutation.

    Attributes:
        query (str): The preprocessed query for fitness evaluation.
        pop_size (int): Total number of individuals in the population.
        pop_selection (int): Number of individuals selected as parents.
        mutation_ratio (float): Probability of mutation for a given individual.
        cutoff (float): Fitness cutoff to record hits.
        n_iterations (int): Number of generations to run.
        maximize (bool): Determines if the algorithm is maximizing or minimizing the fitness.
        sorting (bool): Derived attribute; sorts in ascending order if maximizing is False.
        seed (int): Random seed used for reproducibility.
        bb_manager (BuildingBlockManager): Manager for building block operations.
        fitness_function: Fitness function operator.
        selection_method: Selection method operator.
        crossover_function: Crossover method operator.
        population (List[str]): Current list of peptide sequences.
        run_id (str): Identifier for the current run.
        results_handler (ResultsHandler): Handler for logging and finalizing results.
    """

    def __init__(self,
                 query: str,
                 query_format: str = 'smiles',
                 pop_size: int = 50,
                 pop_selection: int = 10,
                 mutation_ratio: float = 0.8,
                 cutoff: float = 0.5,
                 fitness_function: str = 'map4c',
                 selection_strategy: str = 'greedy',
                 crossover_method: str = 'single_point',
                 n_iterations: int = 1000,
                 run_id: str = 'run',
                 maximize: bool = True,
                 seed: int = 0
                 ):
        """
        Initialize the PDGA algorithm.

        Sets up the algorithm parameters, initializes building block management, retrieves modular operators,
        processes the query, creates an initial population, and logs the configuration.

        Args:
            query (str): The target query used for fitness evaluation.
            query_format (str, optional): Format of the query. Defaults to 'smiles'.
            pop_size (int, optional): Population size. Defaults to 50.
            pop_selection (int, optional): Number of individuals selected as parents. Defaults to 10.
            mutation_ratio (float, optional): Ratio of individuals to be mutated. Defaults to 0.5.
            cutoff (float, optional): Fitness cutoff to record hits. Defaults to 0.5.
            fitness_function (str, optional): Identifier for the fitness function to use. Defaults to 'map4c'.
            selection_strategy (str, optional): Identifier for the selection method. Defaults to 'maximize'.
            crossover_method (str, optional): Identifier for the crossover method. Defaults to 'single_point'.
            n_iterations (int, optional): Number of generations to run. Defaults to 1000.
            run_id (str, optional): Identifier for the current run. Defaults to 'run'.
            maximize (bool, optional): Whether the goal is to maximize the fitness function. Defaults to True.
            seed (int, optional): Random seed for reproducibility. Defaults to 0.
        """
        # Set random seed for reproducibility.
        random.seed(seed)

        # Initialize the PDGA configuration.
        self.pop_size = pop_size
        self.pop_selection = pop_selection
        self.mutation_ratio = mutation_ratio
        self.cutoff = cutoff
        self.n_iterations = n_iterations
        self.maximize = maximize
        self.sorting = not maximize
        self.seed = seed

        # Initialize the building block manager and process the query.
        self.bb_manager = BuildingBlockManager()

        # Retrieve the fitness function and its preprocessing functions, selection, and crossover.
        self.fitness_function = get_fitness_function(fitness_function)
        self.selection_method = get_selection_method(selection_strategy, maximize=self.maximize)
        self.crossover_function = get_crossover_method(crossover_method)

        # Preprocess the query and initialize the population.
        self.query = self.fitness_function.process_query(query, query_format)
        self.population: List[str] = [self.bb_manager.random_linear_seq() for _ in range(self.pop_size)]

        # Initialize the results handler with the desired sort order and run ID.
        timestamp = time.strftime('%y%m%d')
        self.run_id = f'{run_id}_{fitness_function}_{timestamp}_{seed}'
        self.results_handler = ResultsHandler(sort_ascending=self.sorting, run_id=self.run_id)

        # Log the configuration and building block information.
        run_params = {
            'query': self.query,
            'query_format': query_format,
            'pop_size': self.pop_size,
            'pop_selection': self.pop_selection,
            'mutation_ratio': self.mutation_ratio,
            'cutoff': self.cutoff,
            'fitness_function': fitness_function,
            'selection_strategy': selection_strategy,
            'crossover_method': crossover_method,
            'n_iterations': self.n_iterations,
            'run_id': self.run_id,
            'maximize': self.maximize,
            'seed': self.seed,
            'bb_list_count': self.bb_manager.bb_list,
            'ncap_list_count': self.bb_manager.ncap_list,
            'branch_list_count': self.bb_manager.branch_list
        }
        self.results_handler.log_config(run_params)

    def evaluate_fitness(self) -> List[float]:
        """
        Evaluate the fitness for each individual in the current population.

        Iterates through the population, processes each individual, and computes its fitness score
        using the selected fitness function.

        Returns:
            List[float]: A list containing the fitness scores corresponding to each individual.
        """
        fitness_scores: List[float] = []
        for individual in self.population:
            processed_individual = self.fitness_function.process(individual)
            score = self.fitness_function.fitness(processed_individual, self.query)
            fitness_scores.append(score)
        return fitness_scores

    def select_parents(self, fitness_scores: List[float]) -> List[str]:
        """
        Select parent individuals based on their fitness scores.

        Uses the selection method operator to choose a subset of the population as parents.

        Args:
            fitness_scores (List[float]): Fitness scores for the current population.

        Returns:
            List[str]: A list of selected parent sequences.
        """
        return self.selection_method(fitness_scores, self.population, self.pop_selection)

    def apply_crossover(self, parents: List[str]) -> List[str]:
        """
        Generate new offspring using the crossover operator on selected parents.

        Iterates through pairs of parent sequences and applies the crossover function to produce offspring.
        If there is an odd number of parents, the last parent is paired with itself.

        Args:
            parents (List[str]): A list of parent sequences.

        Returns:
            List[str]: A list of offspring sequences generated from the parents.
        """
        next_population: List[str] = []
        for i in range(0, len(parents), 2):
            parent1 = parents[i]
            parent2 = parents[i + 1] if i + 1 < len(parents) else parents[i]
            offspring1, offspring2 = self.crossover_function(parent1, parent2)
            next_population.extend([offspring1, offspring2])
        return next_population

    def apply_mutations(self, population: List[str]) -> List[str]:
        """
        Apply mutations to a population based on the mutation ratio.

        For each individual in the input population, a random chance determines whether the mutation
        operator is applied. If so, the individual is mutated using the building block manager.

        Args:
            population (List[str]): List of sequences to potentially mutate.

        Returns:
            List[str]: The resulting list of sequences after mutation.
        """
        mutated_population: List[str] = []
        for individual in population:
            if random.random() < self.mutation_ratio:
                mutated_individual = mutate(individual, self.bb_manager.bb_list, self.bb_manager.ncap_list, self.bb_manager.branch_list)
                mutated_population.append(mutated_individual)
            else:
                mutated_population.append(individual)
        return mutated_population

    def store_hits(self, fitness_scores: List[float], generation: int) -> None:
        """
        Evaluate and store hits from the current population based on the fitness cutoff.

        For maximization problems, stores an individual if its fitness score is greater than or equal to the cutoff.
        For minimization problems, stores an individual if its fitness score is less than or equal to the cutoff.
        The SMILES representation of each hit is generated and the hit is logged via the ResultsHandler.

        Args:
            fitness_scores (List[float]): Fitness scores for the current population.
            generation (int): The current generation number.

        Returns:
            None
        """
        if self.maximize:
            for individual, score in zip(self.population, fitness_scores):
                if score >= self.cutoff:
                    hit_smiles = self.bb_manager.seq_to_smiles(individual)
                    self.results_handler.add_hit(individual, hit_smiles, score, generation)
        else:
            for individual, score in zip(self.population, fitness_scores):
                if score <= self.cutoff:
                    hit_smiles = self.bb_manager.seq_to_smiles(individual)
                    self.results_handler.add_hit(individual, hit_smiles, score, generation)

    def optimize(self) -> None:
        """
        Run the genetic algorithm optimization over a fixed number of generations.

        For each generation, the algorithm performs the following steps:
            1. Evaluates fitness for the current population.
            2. Selects parent individuals based on fitness scores.
            3. Generates offspring using the crossover operator.
            4. Applies mutation to the offspring.
            5. Stores any hits meeting the fitness cutoff.
            6. Fills the population with random sequences if needed.
        After completing all generations, the results are finalized and saved to a CSV file.

        Returns:
            None
        """
        for gen in tqdm(range(self.n_iterations), desc='Generation'):
            fitness_scores = self.evaluate_fitness()
            parents = self.select_parents(fitness_scores)
            offspring = self.apply_mutations(self.apply_crossover(parents))
            
            offspring_fitness = []
            for child in offspring:
                processed = self.fitness_function.process(child)
                score = self.fitness_function.fitness(processed, self.query)
                offspring_fitness.append(score)

            combined = list(zip(self.population, fitness_scores)) + list(zip(offspring, offspring_fitness))
            combined_sorted = sorted(combined, key=lambda x: x[1], reverse=self.maximize)
            self.population = [seq for seq, _ in combined_sorted[:self.pop_size]]
            self.store_hits([score for _, score in combined_sorted[:self.pop_size]], generation=gen)
        
        self.results_handler.finalize_results()