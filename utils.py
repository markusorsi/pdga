import pandas as pd

class ResultsHandler:
    """
    Handles the collection, processing, and export of PDGA result hits.
    
    This class provides methods to add new hit results (each consisting of a sequence, 
    its corresponding SMILES representation, and a fitness score), retrieve the results
    as a sorted pandas DataFrame, and save them to a CSV file.
    
    Attributes:
        hits (List[Tuple[str, str, float]]): A list of result tuples.
        sort_ascending (bool): Whether to sort the results in ascending order (True) or descending (False).
    """
    
    def __init__(self, sort_ascending: bool = False):
        """
        Initialize the ResultsHandler with an empty list of hits and a sort order.
        
        Parameters:
            sort_ascending (bool): If True, fitness scores are sorted in ascending order;
                                   if False, in descending order. Defaults to False.
        """
        self.hits = []
        self.sort_ascending = sort_ascending
    
    def add_hit(self, sequence: str, smiles: str, fitness: float) -> None:
        """
        Add a new hit result to the collection.
        
        Parameters:
            sequence (str): The sequence of the hit.
            smiles (str): The SMILES representation of the sequence.
            fitness (float): The fitness score associated with the hit.
        """
        self.hits.append((sequence, smiles, fitness))
    
    def get_results_dataframe(self) -> pd.DataFrame:
        """
        Retrieve the stored hit results as a pandas DataFrame sorted by fitness.
        
        Returns:
            pd.DataFrame: A DataFrame with columns 'sequence', 'smiles', and 'fitness', sorted 
                          according to the sort_ascending attribute.
        """
        if not self.hits:
            return pd.DataFrame(columns=['sequence', 'smiles', 'fitness'])
        
        df = pd.DataFrame(self.hits, columns=['sequence', 'smiles', 'fitness'])
        df.sort_values(by='fitness', ascending=self.sort_ascending, inplace=True)
        return df
    
    def save_to_csv(self, filename: str) -> None:
        """
        Save the hit results to a CSV file.
        
        Parameters:
            filename (str): The name of the CSV file to which the results should be saved.
        """
        df = self.get_results_dataframe()
        df.to_csv(filename, index=False)
        print(f'Results saved to {filename}')
