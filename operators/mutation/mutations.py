from abc import ABC, abstractmethod
import re
from random import choice
from typing import List

# -------------------------
# Base class for Mutation Operators
# -------------------------

class MutationOperator(ABC):
    
    @abstractmethod
    def applicable(self, seq: str) -> bool:
        """Return True if this operator can be applied to the given sequence."""
        raise NotImplementedError

    @abstractmethod
    def mutate(self, seq: str) -> str:
        """Return the mutated sequence."""
        raise NotImplementedError

# -------------------------
# Building Block (BB) Operators
# -------------------------

class BBDeletion(MutationOperator):
    def applicable(self, seq: str) -> bool:
        # At least two building blocks present.
        return bool(len(re.findall(r'BB[0-9]{3}', seq)) >= 2)

    def mutate(self, seq: str,) -> str:
        matches = list(re.finditer(r'BB[0-9]{3}-', seq))
        replace = choice(matches)
        return seq[:replace.start()] + seq[replace.end():]

class BBInsertion(MutationOperator):
    def __init__(self, bb_list: List[str]):
        self.bb_list = bb_list

    def applicable(self, seq: str) -> bool:
        # Always applicable.
        return bool(seq)

    def mutate(self, seq: str) -> str:
        matches = list(re.finditer(r'-', seq))
        replace = choice(matches)
        return seq[:replace.start()] + f'-{choice(self.bb_list)}-' + seq[replace.end():]

class BBReplacement(MutationOperator):
    def __init__(self, bb_list: List[str]):
        self.bb_list = bb_list

    def applicable(self, seq: str) -> bool:
        # At least one building block present.
        return bool(re.search(r'BB[0-9]{3}', seq))

    def mutate(self, seq: str):
        matches = list(re.finditer(r'BB[0-9]{3}', seq))
        replace = choice(matches)
        return seq[:replace.start()] + choice(self.bb_list) + seq[replace.end():]

# -------------------------
# N-cap Operators
# -------------------------

class NcapDeletion(MutationOperator):
    def applicable(self, seq: str) -> bool:
        # N-cap present.
        return bool(re.search(r'T[0-9]{3}-', seq))

    def mutate(self, seq: str) -> str:
        matches = list(re.finditer(r'T[0-9]{3}-', seq))
        replace = choice(matches)
        return seq[:replace.start()] + seq[replace.end():]

class NcapInsertion(MutationOperator):
    def __init__(self, ncap_list: List[str]):
        self.ncap_list = ncap_list

    def applicable(self, seq: str) -> bool:
        # No cyclization present.
        # No N-cap present.
        has_c = bool(re.search(r'c', seq))
        has_ncap = bool(re.search(r'T', seq))
        return (not has_c) and (not has_ncap)

    def mutate(self, seq: str) -> str:
        return f'{choice(self.ncap_list)}-' + seq

class NcapReplacement(MutationOperator):
    def __init__(self, ncap_list: List[str]):
        self.ncap_list = ncap_list

    def applicable(self, seq: str) -> bool:
        # N-cap present.
        return bool(re.search(r'T[0-9]{3}', seq))

    def mutate(self, seq: str) -> str:
        matches = list(re.finditer(r'T[0-9]{3}', seq))
        replace = choice(matches)
        return seq[:replace.start()] + choice(self.ncap_list) + seq[replace.end():]

# -------------------------
# Branch Operators
# -------------------------

class BranchDeletion(MutationOperator):
    def applicable(self, seq: str) -> bool:
        # Branching point present.
        return bool(re.search(r'b[0-9]{3}-', seq))

    def mutate(self, seq: str) -> str:
        matches = list(re.finditer(r'b[0-9]{3}-', seq))
        replace = choice(matches)
        return seq[:replace.start()] + seq[replace.end():]

class BranchInsertion(MutationOperator):
    def __init__(self, branch_list: List[str]):
        self.branch_list = branch_list

    def applicable(self, seq: str) -> bool:
        # No cyclization present.
        # No branching point present.
        has_c = bool(re.search(r'c', seq))
        has_branch = bool(re.search(r'b', seq))
        return (not has_c) and (not has_branch)

    def mutate(self, seq: str) -> str:
        matches = list(re.finditer(r'-', seq))
        replace = choice(matches)
        return seq[:replace.start()] + f'-{choice(self.branch_list)}-' + seq[replace.end():]

class BranchReplacement(MutationOperator):
    def __init__(self, branch_list: List[str]):
        self.branch_list = branch_list

    def applicable(self, seq: str) -> bool:
        # Branching point present.
        return bool(re.search(r'b[0-9]{3}', seq))

    def mutate(self, seq: str) -> str:
        matches = list(re.finditer(r'b[0-9]{3}', seq))
        replace = choice(matches)
        return seq[:replace.start()] + choice(self.branch_list) + seq[replace.end():]

# -------------------------
# Cyclization Operators (c-to-n)
# -------------------------

class BreakCToN(MutationOperator):
    def applicable(self, seq: str) -> bool:
        # Cyclization present.
        return bool(re.search(r'c', seq))

    def mutate(self, seq: str) -> str:
        return re.sub(r'c-', '', seq)

class CyclizationCToN(MutationOperator):
    # No cyclization present.
    # No branching point present.
    # No N-cap present.
    def applicable(self, seq: str) -> bool:
        has_c = bool(re.search(r'c', seq))
        has_branch = bool(re.search(r'b', seq))
        has_ncap = bool(re.search(r'T', seq))
        return (not has_c) and (not has_branch) and (not has_ncap)

    def mutate(self, seq: str) -> str:
        return f'c-{seq}'

# -------------------------
# Disulfide Operators
# -------------------------

class BreakDisulfide(MutationOperator):
    def applicable(self, seq: str) -> bool:
        # Disulfide bond present.
        return bool(re.search(r's', seq))

    def mutate(self, seq: str) -> str:
        return re.sub(r'(^s-)|(-s$)|(-s)', '', seq)

class CyclizationDisulfide(MutationOperator):
    def applicable(self, seq: str) -> bool:
        # No disulfide bond present.
        return bool(not re.search(r's', seq))

    def mutate(self, seq: str) -> str:
        positions = list(re.finditer(r'-', seq))
        if len(positions) >= 2:
            replace1 = choice(positions)
            seq = seq[:replace1.start()] + f'-s-' + seq[replace1.end():]
            positions = list(re.finditer(r'-', seq))
            replace2 = choice(positions)
            return seq[:replace2.start()] + f'-s-' + seq[replace2.end():]
        return seq