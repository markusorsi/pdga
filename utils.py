import os
import pandas as pd

class ResultsHandler:
    """
    Handles the collection, processing, and export of PDGA result hits.
    
    This class writes each new hit (a tuple of sequence, SMILES, fitness, and generation)
    to a temporary CSV file as they are generated. This avoids keeping all hits in memory 
    during long runs. Once the run is complete, the temporary file is read, sorted 
    (by fitness, in either ascending or descending order), and the final results are saved 
    to a specified CSV file.
    
    Attributes:
        sort_ascending (bool): Whether to sort results in ascending order (True) or descending order (False).
        temp_filename (str): Name of the temporary CSV file used during the run.
    """
    
    def __init__(self, sort_ascending: bool = False, run_id: str = 'run'):
        """
        Initialize the ResultsHandler with a sort order and temporary results file.
        
        Parameters:
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
        Write the contents of the configuration dictionary to a text file in a 'logs' directory.
        The file is named using the run id and has a '.txt' extension.
        
        Parameters:
            config (dict): The configuration dictionary to log.
        """
        with open(self.log_filename, 'w') as f:
            for key, value in run_params.items():
                f.write(f'{key}: {value}\n')

    def add_hit(self, sequence: str, smiles: str, fitness: float, generation: int) -> None:
        """
        Append a new hit to the temporary results file.
        
        Parameters:
            sequence (str): The peptide sequence.
            smiles (str): The SMILES representation of the sequence.
            fitness (float): The fitness score.
            generation (int): The generation number when the hit was recorded.
        """
        with open(self.temp_filename, 'a') as f:
            f.write(f'{sequence},{smiles},{fitness},{generation}\n')

    def finalize_results(self) -> None:
        """
        Read the temporary results file, remove duplicate hits based on sequence (keeping the best fitness per sequence),
        sort the hits by fitness, and save the final CSV file. Then, replace the temporary file with the final file.
        
        Parameters:
            final_filename (str): The filename for the final sorted results CSV.
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