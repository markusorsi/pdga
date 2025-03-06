import re
from random import choice
from typing import List

from .mutations import (
    MutationOperator,
    BBDeletion,
    BBInsertion,
    BBReplacement,
    NcapDeletion,
    NcapInsertion,
    NcapReplacement,
    BranchDeletion,
    BranchInsertion,
    BranchReplacement,
    BreakCToN,
    CyclizationCToN,
    BreakDisulfide,
    CyclizationDisulfide,
)

def mutate(seq: str, bb_list: List[str], ncap_list: List[str], branch_list: List[str]) -> str:
    """
    Build a pool of mutation operators and apply one at random
    if it is applicable to the given sequence.
    """
    operators: List[MutationOperator] = [
        BBDeletion(),
        BBInsertion(bb_list),
        BBReplacement(bb_list),
        NcapDeletion(),
        NcapInsertion(ncap_list),
        NcapReplacement(ncap_list),
        BranchDeletion(),
        BranchInsertion(branch_list),
        BranchReplacement(branch_list),
        BreakCToN(),
        CyclizationCToN(),
        BreakDisulfide(),
        CyclizationDisulfide(),
    ]
    
    # Filter operators by their applicability to the current sequence.
    applicable_ops = [op for op in operators if op.applicable(seq)]
    
    # If there is exactly one BB token, force BBInsertion.
    if len(re.findall(r'BB[0-9]{3}', seq)) == 1:
        applicable_ops = [op for op in operators if isinstance(op, BBInsertion)]
    
    if applicable_ops:
        op = choice(applicable_ops)
        mutated_seq = op.mutate(seq)
        return mutated_seq
    else:
        print("No applicable mutation operator found.")
        return seq
