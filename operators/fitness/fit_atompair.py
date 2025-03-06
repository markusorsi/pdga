from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
from bblocks.bbmanager import BuildingBlockManager
from .fitness import FitnessOperator

bbmanager = BuildingBlockManager()
ap_generator = rdFingerprintGenerator.GetAtomPairGenerator()

class Atompair(FitnessOperator):
    name = "atompair"

    def process_query(self, query: str):
        """
        Process the query to convert it to a fingerprint.

        Parameters:
            query (str): The query SMILES string.

        Returns:
            The fingerprint of the query molecule.
        """
        mol = Chem.MolFromSmiles(query)
        fingerprint = ap_generator.GetFingerprint(mol)
        return fingerprint

    def process(self, sequence: str):
        """
        Process the peptide sequence to convert it to a fingerprint.

        Parameters:
            sequence (str): The peptide sequence to process.
            translation_dict (dict): A dictionary used to translate the sequence to a SMILES representation.

        Returns:
            The fingerprint of the processed molecule.
        """
        smiles = bbmanager.seq_to_smiles(sequence)
        mol = Chem.MolFromSmiles(smiles)
        fingerprint = ap_generator.GetFingerprint(mol)
        return fingerprint

    def fitness(self, individual, query):
        """
        Calculate the fitness score based on the Tanimoto similarity between
        the fingerprint of an individual and the query fingerprint.

        Parameters:
            individual: The fingerprint of the processed individual.
            query: The fingerprint of the processed query.

        Returns:
            float: The Tanimoto similarity score.
        """
        return TanimotoSimilarity(individual, query)
