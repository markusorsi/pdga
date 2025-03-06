import argparse
import importlib.util
import logging
from pdga import PDGA

def load_config(file_path: str) -> dict:
    """
    Dynamically load a configuration file as a Python module and return the CONFIG dictionary.
    
    Parameters:
        file_path (str): Path to the config file.
    
    Returns:
        dict: The CONFIG dictionary defined in the file.
    """
    spec = importlib.util.spec_from_file_location("config", file_path)
    config_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config_module)
    return config_module.CONFIG

def parse_args():
    # First, create a minimal parser to get the config file path.
    initial_parser = argparse.ArgumentParser(add_help=False)
    initial_parser.add_argument("--config", type=str, default="config.py",
                                help="Path to configuration file")
    initial_args, remaining_args = initial_parser.parse_known_args()
    
    # Load the configuration from the specified file.
    CONFIG = load_config(initial_args.config)
    
    # Now, create the full parser using defaults from the loaded CONFIG.
    parser = argparse.ArgumentParser(
        description="Run the Peptide Design Genetic Algorithm (PDGA)"
    )
    parser.add_argument("--config", type=str, default="config.py",
                        help="Path to configuration file")
    parser.add_argument("--query", type=str, default=CONFIG.get("query", "default_query"),
                        help="Target query for fitness evaluation")
    parser.add_argument("--query_format", type=str, default=CONFIG.get("query_format", "smiles"),
                        help="Format of the query (e.g., 'smiles')")
    parser.add_argument("--pop_size", type=int, default=CONFIG.get("pop_size", 50),
                        help="Population size")
    parser.add_argument("--pop_selection", type=int, default=CONFIG.get("pop_selection", 10),
                        help="Number of individuals selected as parents")
    parser.add_argument("--mutation_ratio", type=float, default=CONFIG.get("mutation_ratio", 0.5),
                        help="Ratio of individuals to be mutated")
    parser.add_argument("--cutoff", type=float, default=CONFIG.get("cutoff", 0.5),
                        help="Fitness cutoff to record hits")
    parser.add_argument("--fitness_function", type=str, default=CONFIG.get("fitness_function", "map4c"),
                        help="Identifier for the fitness function to use")
    parser.add_argument("--selection_strategy", type=str, default=CONFIG.get("selection_strategy", "maximize"),
                        help="Identifier for the selection strategy")
    parser.add_argument("--crossover_method", type=str, default=CONFIG.get("crossover_method", "single_point"),
                        help="Identifier for the crossover method")
    parser.add_argument("--n_iterations", type=int, default=CONFIG.get("n_iterations", 1000),
                        help="Number of generations to run")
    parser.add_argument("--run_id", type=str, default=CONFIG.get("run_id", "pdga_run"),
                        help="Identifier for the run (used in output filenames)")
    parser.add_argument("--seed", type=int, default=CONFIG.get("seed", 0),
                        help="Random seed for reproducibility")
    parser.add_argument("--sort_ascending", action="store_true",
                        help="Sort results in ascending order (default is descending)")
    return parser.parse_args()

def main():
    args = parse_args()
    
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    
    pdga_instance = PDGA(
        query=args.query,
        query_format=args.query_format,
        pop_size=args.pop_size,
        pop_selection=args.pop_selection,
        mutation_ratio=args.mutation_ratio,
        cutoff=args.cutoff,
        fitness_function=args.fitness_function,
        selection_strategy=args.selection_strategy,
        crossover_method=args.crossover_method,
        n_iterations=args.n_iterations,
        run_id=args.run_id,
        seed=args.seed,
        sort_ascending=args.sort_ascending
    )
    
    logging.info("Starting optimization...")
    pdga_instance.optimize()
    logging.info("Optimization complete.")

if __name__ == "__main__":
    main()
