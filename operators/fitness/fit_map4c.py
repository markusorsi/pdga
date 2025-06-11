from rdkit import Chem
from mapchiral import mapchiral
from bblocks.bbmanager import BuildingBlockManager
from .fitness import FitnessOperator

# Initialize the building block manager.
bbmanager = BuildingBlockManager()

class MAP4C(FitnessOperator):
    """
    Fitness operator that calculates molecular fingerprints using the MAP4C algorithm.

    This operator processes both queries and peptide sequences by converting them into molecular fingerprints.
    It then computes the fitness score based on the Jaccard similarity between the fingerprints of the individual
    and the query.

    Attributes:
        name (str): Unique identifier for this fitness operator.
    """
    name = "map4c"

    def process_query(self, query: str, query_format: str):
        """
        Process the query to convert it into a molecular fingerprint using MAP4C.

        Depending on the query format, this method converts the query into a SMILES string:
            - 'smiles': Uses the query as-is.
            - 'sequence': Converts a peptide sequence into a SMILES string.
            - 'pdga_sequence': Converts a PDGA sequence into a SMILES string using the building block manager.
        The resulting SMILES is converted into an RDKit molecule, and then a fingerprint is computed using MAP4C.

        Args:
            query (str): The query input, which can be a SMILES string, a peptide sequence, or a PDGA sequence.
            query_format (str): The format of the query (e.g., 'smiles', 'sequence', or 'pdga_sequence').

        Returns:
            The molecular fingerprint of the query.

        Raises:
            ValueError: If an invalid query format is provided.
        """
        if query_format == 'smiles':
            _query = query
        elif query_format == 'sequence':
            _query = Chem.MolToSmiles(Chem.MolFromSequence(query, flavor=1))
        elif query_format == 'pdga_sequence':
            _query = bbmanager.seq_to_smiles(query)
        else:
            raise ValueError(f'Invalid query format: {query_format}')

        mol = Chem.MolFromSmiles(_query)
        fingerprint = mapchiral.encode(mol)
        return fingerprint

    def process(self, sequence: str):
        """
        Process a peptide sequence to convert it into a molecular fingerprint using MAP4C.

        The method translates the peptide sequence into a SMILES string using the building block manager,
        converts it to an RDKit molecule, and computes its fingerprint using MAP4C.

        Args:
            sequence (str): The peptide sequence to process.

        Returns:
            The molecular fingerprint of the processed sequence.
        """
        smiles = bbmanager.seq_to_smiles(sequence)
        mol = Chem.MolFromSmiles(smiles)
        fingerprint = mapchiral.encode(mol)
        return fingerprint

    def fitness(self, individual, query):
        """
        Calculate the fitness score based on the similarity between fingerprints.

        The fitness score is computed using the Jaccard similarity between the fingerprint of the individual
        and the fingerprint of the query. Higher similarity scores indicate a better match.

        Args:
            individual: The fingerprint of the processed individual.
            query: The fingerprint of the processed query.

        Returns:
            float: The Jaccard similarity score between the individual's fingerprint and the query fingerprint.
        """
        return 1 - mapchiral.jaccard_similarity(individual, query)
