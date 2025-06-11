from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
from bblocks.bbmanager import BuildingBlockManager
from .fitness import FitnessOperator

# Initialize the building block manager and atom pair fingerprint generator.
bbmanager = BuildingBlockManager()
ap_generator = rdFingerprintGenerator.GetAtomPairGenerator()

class Atompair(FitnessOperator):
    """
    Fitness operator that calculates molecular fingerprints using the atom pair method.

    This operator processes queries and peptide sequences by converting them into atom pair fingerprints.
    It then computes the fitness score as the Tanimoto similarity between the fingerprints of the individual
    and the query molecule.

    Attributes:
        name (str): Unique identifier for this fitness operator.
    """
    name = "atompair"

    def process_query(self, query: str, query_format: str):
        """
        Process the query to convert it into an atom pair fingerprint.

        Depending on the query format, this method converts the query into a SMILES string:
          - 'smiles': Uses the query directly.
          - 'sequence': Converts a peptide sequence into a SMILES string.
          - 'pdga_sequence': Converts a PDGA sequence into a SMILES string using the building block manager.
        The SMILES is then converted into an RDKit molecule, and its atom pair fingerprint is computed.

        Args:
            query (str): The query input, which may be a SMILES string, a peptide sequence, or a PDGA sequence.
            query_format (str): The format of the query (e.g., 'smiles', 'sequence', or 'pdga_sequence').

        Returns:
            The atom pair fingerprint of the query molecule.

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
        fingerprint = ap_generator.GetFingerprint(mol)
        return fingerprint

    def process(self, sequence: str):
        """
        Process a peptide sequence to convert it into an atom pair fingerprint.

        The method translates the peptide sequence into a SMILES string using the building block manager,
        converts it into an RDKit molecule, and computes its atom pair fingerprint.

        Args:
            sequence (str): The peptide sequence to process.

        Returns:
            The atom pair fingerprint of the processed molecule.
        """
        smiles = bbmanager.seq_to_smiles(sequence)
        mol = Chem.MolFromSmiles(smiles)
        fingerprint = ap_generator.GetFingerprint(mol)
        return fingerprint

    def fitness(self, individual, query):
        """
        Calculate the fitness score based on the Tanimoto similarity between fingerprints.

        The fitness score is computed as the Tanimoto similarity between the fingerprint of the processed
        individual and the fingerprint of the processed query. Higher scores indicate greater similarity.

        Args:
            individual: The atom pair fingerprint of the processed individual.
            query: The atom pair fingerprint of the processed query.

        Returns:
            float: The Tanimoto similarity score between the individual's fingerprint and the query fingerprint.
        """
        return 1 - TanimotoSimilarity(individual, query)
