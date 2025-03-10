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
    Apply a random applicable mutation operator to the given sequence.

    This function builds a pool of mutation operators for building block, N-cap, branch, cyclization,
    and disulfide mutations. It filters the operators by checking which are applicable to the input sequence.
    If exactly one building block is present in the sequence, the operator pool is restricted to only include
    the BBInsertion operator. A random operator from the applicable pool is then selected to mutate the sequence.

    Args:
        seq (str): The input sequence to be mutated.
        bb_list (List[str]): List of available building block identifiers.
        ncap_list (List[str]): List of available N-cap identifiers.
        branch_list (List[str]): List of available branch identifiers.

    Returns:
        str: The mutated sequence if a mutation operator is applicable; otherwise, returns the original sequence.
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
