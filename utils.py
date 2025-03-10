import os
import pandas as pd

class ResultsHandler:
    """
    Handles the collection, processing, and export of PDGA result hits.

    This class writes each new hit (a tuple of sequence, SMILES, fitness, and generation)
    to a temporary CSV file as they are generated. This avoids keeping all hits in memory 
    during long runs. Once the run is complete, the temporary file is read, sorted 
    (by fitness, in either ascending or descending order), duplicates are removed,
    and the final results are saved to a specified CSV file.

    Attributes:
        sort_ascending (bool): Whether to sort results in ascending order (True) or descending order (False).
        logs_dir (str): Directory path for storing log files.
        results_dir (str): Directory path for storing result files.
        log_filename (str): Filename for storing run configuration logs.
        temp_filename (str): Filename for the temporary CSV file used during the run.
        final_filename (str): Filename for the final sorted results CSV.
    """
    
    def __init__(self, sort_ascending: bool = False, run_id: str = 'run'):
        """
        Initializes the ResultsHandler with the specified sort order and run identifier.

        This method sets up directories for logs and results, creates filenames based on the 
        provided run_id, and initializes a temporary results file with CSV headers.

        Args:
            sort_ascending (bool): If True, sort results in ascending order; if False, descending. Defaults to False.
            run_id (str): Identifier for the run, used to create temporary and final filenames.
        """
        # Create directory to store logs.
        self.logs_dir = 'logs'
        if not os.path.exists(self.logs_dir):
            os.makedirs(self.logs_dir)
        
        # Create directory to store results.
        self.results_dir = 'results'
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

        self.sort_ascending = sort_ascending
        self.log_filename = f'{self.logs_dir}/{run_id}_log.txt'
        self.temp_filename = f'{self.results_dir}/{run_id}_temp.txt'
        self.final_filename = f'{self.results_dir}/{run_id}.txt'
        
        # Clear any existing temporary results file by creating a new one with headers.
        with open(self.temp_filename, 'w') as f:
            f.write('sequence,smiles,fitness,generation\n')

    def log_config(self, run_params: dict) -> None:
        """
        Logs the configuration parameters for the run.

        Writes each key-value pair from the configuration dictionary to a text file located
        in the logs directory. The file is named using the run identifier and has a '.txt'
        extension.

        Args:
            run_params (dict): A dictionary containing configuration parameters and their values.

        Returns:
            None
        """
        with open(self.log_filename, 'w') as f:
            for key, value in run_params.items():
                f.write(f'{key}: {value}\n')

    def add_hit(self, sequence: str, smiles: str, fitness: float, generation: int) -> None:
        """
        Appends a new result hit to the temporary CSV file.

        Each hit consists of a peptide sequence, its SMILES representation, a fitness score,
        and the generation number when the hit was recorded. The hit is added as a new line 
        in the temporary file.

        Args:
            sequence (str): The peptide sequence.
            smiles (str): The SMILES representation of the sequence.
            fitness (float): The fitness score of the hit.
            generation (int): The generation number during which the hit was recorded.

        Returns:
            None
        """
        with open(self.temp_filename, 'a') as f:
            f.write(f'{sequence},{smiles},{fitness},{generation}\n')

    def finalize_results(self) -> None:
        """
        Finalizes the results by processing the temporary CSV file.

        The method reads the temporary file into a DataFrame, removes duplicate hits based 
        on the sequence, sorts the hits by fitness in the specified order, deletes the
        temporary file, and saves the final results to a CSV file. A message is printed 
        to indicate where the final results have been saved.

        Returns:
            None

        Raises:
            pd.errors.EmptyDataError: If the temporary file is empty, an empty DataFrame is created.
            OSError: If there is an issue removing the temporary file.
        """
        try:
            dtype = {'sequence': str, 'smiles': str, 'fitness': float, 'generation': int}
            df = pd.read_csv(self.temp_filename, dtype=dtype)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=['sequence', 'smiles', 'fitness', 'generation'])
        
        # Sort and remove duplicates based on sequence, keeping the best fitness.
        df.sort_values(by='fitness', ascending=self.sort_ascending, inplace=True)
        df = df.drop_duplicates(subset='sequence', keep='first')
        
        # Remove the temporary file and save final results.
        os.remove(self.temp_filename)
        df.to_csv(self.final_filename, index=False)
        print(f'Results saved to {self.final_filename}')