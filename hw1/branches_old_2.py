# Maximum likelihood branch length calculation
# Eliot BUsh, Aug 2016
# Joe Wirth, Aug 2022
# Cevi Bainton and Kenneth Mitchell, Sept 2022

from branchHelper import *
from UnrootedTree import Utree
memo_score_one_char = dict()

def branchLenLikelihood(leftSub:tuple, rightSub:tuple, seqD:dict, branchLen:float) -> float:
    """
    Helper function calculates probability of unrooted tree given for branch length. 
    Requires left and right subtrees around selected branch.

    INPUTS:
        • leftSub = left subtree connected to branch (tuple, Utree format)
        • rightSub = right subtree connected to branch (tuple, Utree format)
        • seqD = dictionary of sequences at tips
        • branchLen = proposed branch len (float)

    OUTPUTS:
       • the likelihood given the branch len, calculated by ((probchanged or same) * likelihood_right * likelihoodleft * 1/4)
    """
    total_likelihood_of_seq = 1 #
    likelihood_char_changed = jukesCantorProbability(branchLen) 
    likelihood_char_same = 1 - likelihood_char_changed

    seq_len = len(list(seqD.values())[0]) # get the first sequence from dictionary to learn length

    for pos in range(seq_len): # iterate through every position in the sequence
        likelihood_at_pos = 0

        for lchar in CHARS: # iterate through every possible character for the left tree
            likelihood_left_char = scoreOneChar(leftSub, lchar, pos, seqD, memo_score_one_char) # and score the left tree

            for rchar in CHARS: # iterate through every possible character for the right tree
                likelihood_right_char = scoreOneChar(rightSub, rchar, pos, seqD, memo_score_one_char) # and score the right tree

                if lchar == rchar: # if we didn't change our character, multiply by the chance it stayed the same
                    likelihood_at_pos += likelihood_char_same * likelihood_left_char * likelihood_right_char * 1/4 # add to likelihood

                else: # if we did change our character, multiply by the chance it changed to a specific character
                    likelihood_at_pos += likelihood_char_changed * likelihood_left_char * likelihood_right_char * 1/4 * 1/3 # add to likelihood

        total_likelihood_of_seq *= likelihood_at_pos # multiply to total likelihood

    return total_likelihood_of_seq

    

def getRootedSubtree(tree:Utree, node:int, branch:int) -> tuple:
    # TODO: write this function
    pass
    


def optimizeBranches(tree:Utree, seqD:dict, delta:float) -> tuple[Utree,float]:
    # TODO: write this function
    pass


##### RUNNER for testing #####
def test():
    print('test1')
    __main('example1.nwk','example1.fna',1.01,'example1Out.nwk')
    print()
    print('test2')
    __main('example2.nwk','example2.fna',1.01,'example2Out.nwk')


# WRAPPER FUNCTION
def __main(treeFN:str, fastaFN:str, delta:float, outTreeFN:str) -> None:
    """ main:
            Accepts a tree filename, a fasta filename, a value to adjust branch
            lengths, and a output filename as inputs. Optimizes the branch len-
            gths for the provided tree. Prints the statistics of the tree and
            saves it to the specified file. Does not return.
    """
    # constant
    EOL = "\n"
    
    # read data into memory
    utree = loadTree(treeFN)
    seqD = loadFasta(fastaFN)

    # get some preliminary branch lengths
    treeT = preliminaryBranchLengths(utree.toTuple(),seqD)
    utree = Utree.fromTuple(treeT)
    
    # optimize the tree
    utree,likely = optimizeBranches(utree,seqD,delta)

    # print the tree data for the final tree
    print("Final branch lengths:")
    print(utree, EOL)

    print("Final Likelihood:", format(likely,".5e"))

    # write the final tree to file
    writeUtree(utree,outTreeFN)
    
    print("Tree (arbitrarily rooted) with branch lengths written to:",outTreeFN)
    

seqD = loadFasta('example1.fna')
loadedtree = loadTree('example1.nwk')
temptree  = preliminaryBranchLengths(loadedtree.toTuple(),seqD)
ttree = Utree.fromTuple(temptree)
