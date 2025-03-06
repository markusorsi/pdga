from rdkit import Chem
from mapchiral import mapchiral
from bblocks.bbmanager import BuildingBlockManager
from .fitness import FitnessOperator

bbmanager = BuildingBlockManager()

class MAP4C(FitnessOperator):
    name = "map4c"

    def process_query(self, query: str):
        """
        Process the query to convert it to a fingerprint.
        """
        mol = Chem.MolFromSmiles(query)
        fingerprint = mapchiral.encode(mol)
        return fingerprint

    def process(self, sequence: str):
        """
        Process the peptide sequence to convert it to a fingerprint.
        """
        smiles = bbmanager.seq_to_smiles(sequence)
        mol = Chem.MolFromSmiles(smiles)
        fingerprint = mapchiral.encode(mol)
        return fingerprint

    def fitness(self, individual, query):
        """
        Calculate the fitness score based on the similarity to the query.
        """
        return mapchiral.jaccard_similarity(individual, query)
