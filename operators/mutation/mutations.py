from abc import ABC, abstractmethod
import re
from random import choice
from typing import List

# -------------------------
# Base class for Mutation Operators
# -------------------------

class MutationOperator(ABC):
    """
    Abstract base class for mutation operators.

    MutationOperator defines the interface for mutation operators used to modify sequences
    in the PDGA algorithm. Subclasses must implement both the `applicable` and `mutate` methods.
    """

    @abstractmethod
    def applicable(self, seq: str) -> bool:
        """
        Determine if the mutation operator can be applied to the given sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if the operator can be applied; False otherwise.
        """
        raise NotImplementedError

    @abstractmethod
    def mutate(self, seq: str) -> str:
        """
        Mutate the given sequence and return the modified sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The mutated sequence.
        """
        raise NotImplementedError

# -------------------------
# Building Block (BB) Operators
# -------------------------

class BBDeletion(MutationOperator):
    """
    Mutation operator for deleting a building block from a sequence.

    The operator checks if there are at least two building blocks (formatted as 'BBxxx')
    present in the sequence. If so, it randomly removes one occurrence (including a trailing hyphen).
    """

    def applicable(self, seq: str) -> bool:
        """
        Check if the sequence has at least two building blocks.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if at least two building blocks are found; False otherwise.
        """
        return bool(len(re.findall(r'BB[0-9]{3}', seq)) >= 2)

    def mutate(self, seq: str) -> str:
        """
        Delete a randomly selected building block from the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with one building block removed.
        """
        matches = list(re.finditer(r'BB[0-9]{3}-', seq))
        replace = choice(matches)
        return seq[:replace.start()] + seq[replace.end():]

class BBInsertion(MutationOperator):
    """
    Mutation operator for inserting a building block into a sequence.

    This operator inserts a new building block (randomly selected from a provided list)
    into the sequence at a random hyphen position.
    """

    def __init__(self, bb_list: List[str]):
        """
        Initialize BBInsertion with a list of available building blocks.

        Args:
            bb_list (List[str]): List of building block identifiers.
        """
        self.bb_list = bb_list

    def applicable(self, seq: str) -> bool:
        """
        Determine if insertion is applicable to the sequence.

        Since insertion is always possible if the sequence exists, this returns True
        when the sequence is non-empty.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if the sequence is non-empty; False otherwise.
        """
        return bool(seq)

    def mutate(self, seq: str) -> str:
        """
        Insert a random building block into the sequence at a random hyphen position.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with an additional building block inserted.
        """
        matches = list(re.finditer(r'-', seq))
        replace = choice(matches)
        return seq[:replace.start()] + f'-{choice(self.bb_list)}-' + seq[replace.end():]

class BBReplacement(MutationOperator):
    """
    Mutation operator for replacing a building block in a sequence.

    This operator randomly selects a building block occurrence (formatted as 'BBxxx')
    in the sequence and replaces it with a random building block from a provided list.
    """

    def __init__(self, bb_list: List[str]):
        """
        Initialize BBReplacement with a list of available building blocks.

        Args:
            bb_list (List[str]): List of building block identifiers.
        """
        self.bb_list = bb_list

    def applicable(self, seq: str) -> bool:
        """
        Check if there is at least one building block present in the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if a building block is found; False otherwise.
        """
        return bool(re.search(r'BB[0-9]{3}', seq))

    def mutate(self, seq: str) -> str:
        """
        Replace a randomly selected building block with another from the available list.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with one building block replaced.
        """
        matches = list(re.finditer(r'BB[0-9]{3}', seq))
        replace = choice(matches)
        return seq[:replace.start()] + choice(self.bb_list) + seq[replace.end():]

# -------------------------
# N-cap Operators
# -------------------------

class NcapDeletion(MutationOperator):
    """
    Mutation operator for deleting an N-cap from a sequence.

    The operator looks for an N-cap pattern formatted as 'Txxx-' and removes it.
    """

    def applicable(self, seq: str) -> bool:
        """
        Check if an N-cap is present in the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if an N-cap is found; False otherwise.
        """
        return bool(re.search(r'T[0-9]{3}-', seq))

    def mutate(self, seq: str) -> str:
        """
        Delete a randomly selected N-cap from the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with one N-cap removed.
        """
        matches = list(re.finditer(r'T[0-9]{3}-', seq))
        replace = choice(matches)
        return seq[:replace.start()] + seq[replace.end():]

class NcapInsertion(MutationOperator):
    """
    Mutation operator for inserting an N-cap into a sequence.

    This operator inserts a new N-cap (randomly selected from a provided list)
    at the beginning of the sequence if no cyclization or N-cap is already present.
    """

    def __init__(self, ncap_list: List[str]):
        """
        Initialize NcapInsertion with a list of available N-caps.

        Args:
            ncap_list (List[str]): List of N-cap identifiers.
        """
        self.ncap_list = ncap_list

    def applicable(self, seq: str) -> bool:
        """
        Determine if an N-cap can be inserted.

        Insertion is applicable if the sequence does not contain cyclization ('c')
        and does not already have an N-cap ('T').

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if insertion is allowed; False otherwise.
        """
        has_c = bool(re.search(r'c', seq))
        has_ncap = bool(re.search(r'T', seq))
        return (not has_c) and (not has_ncap)

    def mutate(self, seq: str) -> str:
        """
        Insert a random N-cap at the beginning of the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with a new N-cap prepended.
        """
        return f'{choice(self.ncap_list)}-' + seq

class NcapReplacement(MutationOperator):
    """
    Mutation operator for replacing an existing N-cap in a sequence.

    This operator identifies an N-cap pattern (formatted as 'Txxx') in the sequence and
    replaces it with a random N-cap from a provided list.
    """

    def __init__(self, ncap_list: List[str]):
        """
        Initialize NcapReplacement with a list of available N-caps.

        Args:
            ncap_list (List[str]): List of N-cap identifiers.
        """
        self.ncap_list = ncap_list

    def applicable(self, seq: str) -> bool:
        """
        Check if an N-cap is present in the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if an N-cap is found; False otherwise.
        """
        return bool(re.search(r'T[0-9]{3}', seq))

    def mutate(self, seq: str) -> str:
        """
        Replace a randomly selected N-cap in the sequence with another N-cap.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with an N-cap replaced.
        """
        matches = list(re.finditer(r'T[0-9]{3}', seq))
        replace = choice(matches)
        return seq[:replace.start()] + choice(self.ncap_list) + seq[replace.end():]

# -------------------------
# Branch Operators
# -------------------------

class BranchDeletion(MutationOperator):
    """
    Mutation operator for deleting a branch from a sequence.

    The operator detects a branch pattern (formatted as 'bxxx-') and removes it.
    """

    def applicable(self, seq: str) -> bool:
        """
        Check if a branching point is present in the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if a branch is found; False otherwise.
        """
        return bool(re.search(r'b[0-9]{3}-', seq))

    def mutate(self, seq: str) -> str:
        """
        Delete a randomly selected branch from the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with a branch removed.
        """
        matches = list(re.finditer(r'b[0-9]{3}-', seq))
        replace = choice(matches)
        return seq[:replace.start()] + seq[replace.end():]

class BranchInsertion(MutationOperator):
    """
    Mutation operator for inserting a branch into a sequence.

    This operator inserts a new branch (randomly selected from a provided list) into the sequence
    at a random hyphen position, provided that there is no cyclization or existing branch.
    """

    def __init__(self, branch_list: List[str]):
        """
        Initialize BranchInsertion with a list of available branches.

        Args:
            branch_list (List[str]): List of branch identifiers.
        """
        self.branch_list = branch_list

    def applicable(self, seq: str) -> bool:
        """
        Determine if a branch can be inserted.

        Insertion is allowed if the sequence does not contain cyclization ('c') and no branch ('b') exists.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if insertion is allowed; False otherwise.
        """
        has_c = bool(re.search(r'c', seq))
        has_branch = bool(re.search(r'b', seq))
        return (not has_c) and (not has_branch)

    def mutate(self, seq: str) -> str:
        """
        Insert a random branch into the sequence at a random hyphen position.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with a new branch inserted.
        """
        matches = list(re.finditer(r'-', seq))
        replace = choice(matches)
        return seq[:replace.start()] + f'-{choice(self.branch_list)}-' + seq[replace.end():]

class BranchReplacement(MutationOperator):
    """
    Mutation operator for replacing a branch in a sequence.

    This operator randomly selects an existing branch (formatted as 'bxxx') in the sequence
    and replaces it with a branch from a provided list.
    """

    def __init__(self, branch_list: List[str]):
        """
        Initialize BranchReplacement with a list of available branches.

        Args:
            branch_list (List[str]): List of branch identifiers.
        """
        self.branch_list = branch_list

    def applicable(self, seq: str) -> bool:
        """
        Check if a branch is present in the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if a branch is found; False otherwise.
        """
        return bool(re.search(r'b[0-9]{3}', seq))

    def mutate(self, seq: str) -> str:
        """
        Replace a randomly selected branch with another branch.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with a branch replaced.
        """
        matches = list(re.finditer(r'b[0-9]{3}', seq))
        replace = choice(matches)
        return seq[:replace.start()] + choice(self.branch_list) + seq[replace.end():]

# -------------------------
# Cyclization Operators (c-to-n)
# -------------------------

class BreakCToN(MutationOperator):
    """
    Mutation operator for breaking cyclization (c-to-n) in a sequence.

    This operator detects cyclization represented by 'c' and removes the pattern.
    """

    def applicable(self, seq: str) -> bool:
        """
        Check if cyclization is present in the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if a cyclization pattern ('c') is found; False otherwise.
        """
        return bool(re.search(r'c', seq))

    def mutate(self, seq: str) -> str:
        """
        Remove cyclization from the sequence by deleting the pattern 'c-'.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with cyclization removed.
        """
        return re.sub(r'c-', '', seq)

class CyclizationCToN(MutationOperator):
    """
    Mutation operator for introducing cyclization (c-to-n) into a sequence.

    This operator adds cyclization to a sequence that currently has no cyclization ('c'),
    no branch ('b'), and no N-cap ('T').
    """

    def applicable(self, seq: str) -> bool:
        """
        Determine if cyclization can be applied to the sequence.

        Cyclization is applicable if the sequence lacks cyclization, branches, and N-caps.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if cyclization can be introduced; False otherwise.
        """
        has_c = bool(re.search(r'c', seq))
        has_branch = bool(re.search(r'b', seq))
        has_ncap = bool(re.search(r'T', seq))
        return (not has_c) and (not has_branch) and (not has_ncap)

    def mutate(self, seq: str) -> str:
        """
        Introduce cyclization by prepending 'c-' to the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with cyclization introduced.
        """
        return f'c-{seq}'

# -------------------------
# Disulfide Operators
# -------------------------

class BreakDisulfide(MutationOperator):
    """
    Mutation operator for breaking a disulfide bond in a sequence.

    This operator detects the disulfide bond represented by 's' and removes its occurrence.
    """

    def applicable(self, seq: str) -> bool:
        """
        Check if a disulfide bond is present in the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if a disulfide bond ('s') is found; False otherwise.
        """
        return bool(re.search(r's', seq))

    def mutate(self, seq: str) -> str:
        """
        Remove the disulfide bond from the sequence.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with the disulfide bond removed.
        """
        return re.sub(r'(s-)|(-s)', '', seq)

class CyclizationDisulfide(MutationOperator):
    """
    Mutation operator for forming a disulfide bond (cyclization) in a sequence.

    This operator introduces disulfide bonds (represented by 's') into the sequence by inserting
    the bond at two different positions, provided there are at least two hyphens in the sequence.
    """

    def applicable(self, seq: str) -> bool:
        """
        Determine if a disulfide bond can be formed in the sequence.

        The bond can be introduced if no disulfide bond is currently present.

        Args:
            seq (str): The input sequence.

        Returns:
            bool: True if a disulfide bond is not present; False otherwise.
        """
        return bool(not re.search(r's', seq))

    def mutate(self, seq: str) -> str:
        """
        Introduce disulfide bonds into the sequence at two random hyphen positions.

        Args:
            seq (str): The input sequence.

        Returns:
            str: The sequence with disulfide bonds inserted if applicable; otherwise, the original sequence.
        """
        positions = list(re.finditer(r'-', seq))
        if len(positions) >= 2:
            replace1 = choice(positions)
            seq = seq[:replace1.start()] + f'-s-' + seq[replace1.end():]
            positions = list(re.finditer(r'-', seq))
            replace2 = choice(positions)
            return seq[:replace2.start()] + f'-s-' + seq[replace2.end():]
        return seq
