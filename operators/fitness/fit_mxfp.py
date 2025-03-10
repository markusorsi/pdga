from rdkit import Chem
from mxfp import mxfp
from bblocks.bbmanager import BuildingBlockManager
from .fitness import FitnessOperator

# Initialize building block manager and MXFP calculator.
bbmanager = BuildingBlockManager()
mxfp_calculator = mxfp.MXFPCalculator()

class MXFP(FitnessOperator):
    """
    Fitness operator that calculates the MXFP fingerprint.

    This operator processes queries and peptide sequences by converting them into their
    MXFP fingerprints, and then calculates the fitness score as the sum of absolute differences 
    (Manhattan distance) between the fingerprints.

    Attributes:
        name (str): Unique identifier for this fitness operator.
    """
    name = "mxfp"

    def process_query(self, query: str, query_format: str):
        """
        Process the query to convert it into a molecular fingerprint.

        Depending on the query format, this method converts the query into a SMILES string,
        then creates an RDKit molecule and computes its fingerprint using the MXFP calculator.

        Args:
            query (str): The query input, which can be a SMILES string, a sequence, or a PDGA sequence.
            query_format (str): The format of the query. Expected values are 'smiles', 'sequence', or 'pdga_sequence'.

        Returns:
            The fingerprint of the query molecule.

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
        fingerprint = mxfp_calculator.mxfp_from_mol(mol)
        return fingerprint

    def process(self, sequence: str):
        """
        Process a peptide sequence by converting it into a molecular fingerprint using MXFP.

        The method translates the peptide sequence into a SMILES string using the building block
        manager, then generates an RDKit molecule and computes its fingerprint.

        Args:
            sequence (str): The peptide sequence to process.

        Returns:
            The fingerprint of the processed molecule.
        """
        smiles = bbmanager.seq_to_smiles(sequence)
        mol = Chem.MolFromSmiles(smiles)
        fingerprint = mxfp_calculator.mxfp_from_mol(mol)
        return fingerprint

    def fitness(self, individual, query):
        """
        Calculate the fitness score based on the difference between the individual's fingerprint and the query fingerprint.

        The fitness score is computed as the sum of absolute differences between corresponding
        elements of the individual's fingerprint and the query fingerprint (Manhattan distance).
        Lower scores indicate higher similarity between the fingerprints.

        Args:
            individual: The fingerprint of the processed individual.
            query: The fingerprint of the processed query.

        Returns:
            float: The similarity score computed as the sum of absolute differences.
        """
        return sum(abs(individual - query))
