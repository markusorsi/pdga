from rdkit import Chem
from mxfp import mxfp
from bblocks.bbmanager import BuildingBlockManager
from .fitness import FitnessOperator

bbmanager = BuildingBlockManager()
mxfp_calculator = mxfp.MXFPCalculator()

class MXFP(FitnessOperator):
    name = "mxfp"

    def process_query(self, query: str):
        """
        Process the query to convert it to a fingerprint using MXFP.

        Parameters:
            query (str): The query SMILES string.

        Returns:
            The fingerprint of the query molecule.
        """
        mol = Chem.MolFromSmiles(query)
        fingerprint = mxfp_calculator.mxfp_from_mol(mol)
        return fingerprint

    def process(self, sequence: str):
        """
        Process the peptide sequence to convert it to a fingerprint using MXFP.

        Parameters:
            sequence (str): The peptide sequence to process.
            translation_dict (dict): A dictionary to translate the sequence into a SMILES string.

        Returns:
            The fingerprint of the processed molecule.
        """
        smiles = bbmanager.seq_to_smiles(sequence)
        mol = Chem.MolFromSmiles(smiles)
        fingerprint = mxfp_calculator.mxfp_from_mol(mol)
        return fingerprint

    def fitness(self, individual, query):
        """
        Calculate the fitness score based on the difference between the individual’s fingerprint and the query fingerprint.
        
        Parameters:
            individual: The fingerprint of the processed individual.
            query: The fingerprint of the processed query.
        
        Returns:
            float: The similarity score computed as the sum of absolute differences.
        """
        return sum(abs(individual - query))
