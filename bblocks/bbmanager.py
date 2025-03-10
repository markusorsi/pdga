from random import choice, randint
import pandas as pd

class BuildingBlockManager:
    """
    Manages building block data for the PDGA algorithm.

    This class loads building block data from CSV files, creates lists of available building block IDs,
    and merges translation dictionaries for converting sequence tokens into their corresponding SMILES fragments.
    It also provides utility methods to generate random sequences and to translate a sequence into a SMILES string.

    Attributes:
        bb_list (List[str]): List of building block IDs from the main building blocks CSV.
        ncap_list (List[str]): List of N-cap IDs.
        branch_list (List[str]): List of branch IDs.
        translation_dict (dict): Dictionary mapping building block IDs (and custom tokens 'c' and 's') to SMILES strings.
    """

    def __init__(self,
                 monomers_csv: str = 'bblocks/bb_monomers.csv',
                 ncaps_csv: str = 'bblocks/bb_ncaps.csv',
                 branches_csv: str = 'bblocks/bb_branches.csv',
                 additional_csv: str = 'bblocks/bb_additional.csv'):
        """
        Initialize the BuildingBlockManager by loading CSV files containing building block data.

        Reads in four CSV files containing different sets of building block information and constructs
        translation dictionaries for each. The individual dictionaries are then merged into a single dictionary.
        Custom translations for special tokens are added afterwards.

        Args:
            monomers_csv (str): Path to the CSV file containing monomer building blocks.
            ncaps_csv (str): Path to the CSV file containing N-cap building blocks.
            branches_csv (str): Path to the CSV file containing branch building blocks.
            additional_csv (str): Path to the CSV file containing additional building block data.
        """
        # Load CSV files.
        building_blocks = pd.read_csv(monomers_csv)
        ncaps = pd.read_csv(ncaps_csv)
        branches = pd.read_csv(branches_csv)
        additional = pd.read_csv(additional_csv)

        # Create individual translation dictionaries.
        bb_dict = dict(zip(building_blocks.ID, building_blocks.SMILES))
        ncap_dict = dict(zip(ncaps.ID, ncaps.SMILES))
        branch_dict = dict(zip(branches.ID, branches.SMILES))
        additional_dict = dict(zip(additional.ID, additional.SMILES))

        # Store lists of building block IDs.
        self.bb_list = building_blocks.ID.values.tolist()
        self.ncap_list = ncaps.ID.values.tolist()
        self.branch_list = branches.ID.values.tolist()

        # Merge all translation dictionaries.
        self.translation_dict = {**bb_dict, **ncap_dict, **branch_dict, **additional_dict}

        # Add custom translations for special tokens.
        self.translation_dict['c'] = '9'
        self.translation_dict['s'] = 'NC(CS7)C(=O)'

    def random_linear_seq(self, min_len: int = 5, max_len: int = 20) -> str:
        """
        Generate a random linear sequence from available building blocks.

        The sequence is created by randomly selecting building block IDs from the bb_list and joining
        them with hyphens. The number of building blocks in the sequence is randomly determined
        between the specified minimum and maximum lengths.

        Args:
            min_len (int): Minimum number of building blocks to include in the sequence.
            max_len (int): Maximum number of building blocks to include in the sequence.

        Returns:
            str: A randomly generated sequence of building block IDs separated by hyphens.
        """
        seq = choice(self.bb_list)
        for _ in range(randint(min_len, max_len)):
            seq += '-' + choice(self.bb_list)
        return seq

    def seq_to_smiles(self, seq: str) -> str:
        """
        Translate a sequence of building block IDs into a SMILES string.

        Splits the input sequence on hyphens and converts each token to its corresponding SMILES fragment
        using the translation_dict. The fragments are concatenated to form a complete SMILES string.
        Additional modifications are applied based on the presence of certain tokens:
          - If the sequence contains 'b', appends '8' to the SMILES string to close the branching ring.
          - If the sequence contains 'c', rearranges parts of the SMILES string and appends '9' to close the C-to-N ring.
          - Otherwise, appends 'O' to the SMILES string.

        Args:
            seq (str): A sequence of building block IDs separated by hyphens.

        Returns:
            str: The SMILES string corresponding to the input sequence. Returns an empty string if an error occurs.
        """
        smiles = ''
        try:
            for element in seq.split('-'):
                smiles += self.translation_dict.get(element, '')
            if 'b' in seq:
                return smiles + '8'
            elif 'c' in seq:
                return smiles[1] + smiles[0] + smiles[2:] + '9'
            else:
                return smiles + 'O'
        except Exception as e:
            print(f'Error processing sequence {seq}: {e}')
            return ''
