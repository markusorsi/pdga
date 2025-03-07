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
    Results (hits) are handled using a dedicated ResultsHandler.
    """

    def __init__(self,
                 query: str,
                 query_format: str = 'smiles',
                 pop_size: int = 50,
                 pop_selection: int = 10,
                 mutation_ratio: float = 0.5,
                 cutoff: float = 0.5,
                 fitness_function: str = 'map4c',
                 selection_strategy: str = 'maximize',
                 crossover_method: str = 'single_point',
                 n_iterations: int = 1000,
                 run_id: str = 'run',
                 maximize: bool = True,
                 seed: int = 0
                 ):
        """
        Initialize the PDGA algorithm.

        :param query: The target query used for fitness evaluation.
        :param query_format: Format of the query (default 'smiles').
        :param pop_size: Population size.
        :param pop_selection: Number of individuals selected as parents.
        :param mutation_ratio: Ratio of individuals to be mutated.
        :param cutoff: Fitness cutoff to record hits.
        :param fitness_function: Identifier for the fitness function to use.
        :param selection_strategy: Identifier for the selection method.
        :param crossover_method: Identifier for the crossover method.
        :param n_iterations: Number of generations to run.
        :param run_id: Identifier for the current run.
        :param seed: Random seed for reproducibility.
        :param maximize: Whether the goal is to maximize the fitness function.
        """
        # Set random seed for reproducibility.
        random.seed(seed)

        # Initialize the PDGA configuration.
        self.query = query
        self.query_format = query_format
        self.pop_size = pop_size
        self.pop_selection = pop_selection
        self.mutation_ratio = mutation_ratio
        self.cutoff = cutoff
        self.n_iterations = n_iterations
        self.maximize = maximize
        self.sort_ascending = not maximize
        self.seed = seed

        # Initialize the building block manager and retrieve building block data.
        self.bb_manager = BuildingBlockManager()

        # Retrieve the fitness function and its preprocessing functions, selection, and crossover.
        self.fitness_function = get_fitness_function(fitness_function)
        self.selection_method = get_selection_method(selection_strategy)
        self.crossover_function = get_crossover_method(crossover_method)

        # Preprocess the query and initialize the population.
        self.query_processed = self.fitness_function.process_query(self.query)
        self.population: List[str] = [self.bb_manager.random_linear_seq() for _ in range(self.pop_size)]

        # Initialize the results handler with the desired sort order and run ID.
        timestamp = time.strftime('%y%m%d')
        self.run_id = f'{run_id}_{fitness_function}_{timestamp}_{seed}'
        self.results_handler = ResultsHandler(sort_ascending=self.sort_ascending, run_id=self.run_id)

        # Log the configuration and building block information.
        run_params = {
            'query': self.query,
            'query_format': self.query_format,
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
        Evaluate the fitness for each individual in the population.

        :return: List of fitness scores.
        """
        fitness_scores: List[float] = []
        for individual in self.population:
            processed_individual = self.fitness_function.process(individual)
            score = self.fitness_function.fitness(processed_individual, self.query_processed)
            fitness_scores.append(score)
        return fitness_scores

    def select_parents(self, fitness_scores: List[float]) -> List[str]:
        """
        Select parent individuals based on fitness scores.

        :param fitness_scores: Fitness scores for the current population.
        :return: List of selected parent sequences.
        """
        return self.selection_method(fitness_scores, self.population, self.pop_selection)

    def apply_crossover(self, parents: List[str]) -> List[str]:
        """
        Apply the crossover operator to generate new offspring.

        :param parents: List of parent sequences.
        :return: List of offspring sequences.
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
        Apply mutation to the population based on the mutation ratio.

        :param population: List of sequences to potentially mutate.
        :return: List of sequences after mutation.
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
        Process the current population and store sequences whose fitness meets the cutoff criterion.
        
        For a maximization problem, a hit is stored if the score is greater than or equal to the cutoff.
        For a minimization problem, a hit is stored if the score is less than or equal to the cutoff.
        
        :param fitness_scores: Fitness scores corresponding to the current population.
        :param generation: The current generation number.
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
        Run the genetic algorithm optimization for the specified number of iterations.
        
        At each generation, fitness is evaluated, parents are selected, offspring are generated via crossover,
        mutations are applied, and hits are stored (along with their generation number). After all generations,
        the results are finalized and saved to a CSV file.
        """
        for gen in tqdm(range(self.n_iterations), desc='Generation'):
            fitness_scores = self.evaluate_fitness()
            parents = self.select_parents(fitness_scores)
            offspring = self.apply_crossover(parents)
            new_population = self.apply_mutations(offspring)
            self.store_hits(fitness_scores, generation=gen)
            
            # Fill the rest of the population with random sequences if needed.
            while len(new_population) < self.pop_size:
                new_population.append(self.bb_manager.random_linear_seq())
            self.population = new_population
        
        self.results_handler.finalize_results()
