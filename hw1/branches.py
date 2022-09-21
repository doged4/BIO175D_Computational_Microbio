# Maximum likelihood branch length calculation
# Eliot Bush, Aug 2016
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
       • a float, the likelihood given the branch length, calculated by ((probchanged or same) * likelihood_right * likelihoodleft * 1/4)
    """
    total_likelihood_of_seq = 1 # initial likelihood is 1
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
    """"
    Returns subtree of Utree object rooted at supplied node, including supplied branch

    INPUTS:
        • tree = a Utree object
        • node = the index of a node in the tree which will be the root node
        • branch = the index of a branch in tree which will be the leading branch of the root
    OUTPUTS:
        • a tuple, a rooted subtree of tree where node is the root node and branch is the leading branch
    """
    # Get data from tree for final outpus
    nName = tree.getName(node)
    bLength = tree.getBranchLength(branch)

    # BASE CASE: if node is a leaf
    if tree.isLeaf(node):
        newTree = ((nName, bLength),
                 (),
                 ())
        return newTree

    neighboringBranches = tree.getBranches(node) # Get all branches attached to node of interest including given branch

    newTreeList = [] # This will contain all rooted subtrees we make, which should be two

    # RECURSIVE CASE: recurse through each neighboring branch from rooting node, make new rooted subtrees out of each new branch
    for newBranch in neighboringBranches:
        if newBranch != branch: # The overall rooted tree should start with the argument branch, no subtree should include it
            
            neighboringNodes = tree.getNodes(newBranch)

            for newNode in neighboringNodes: # Get a newNode that is a neighbor of our original node via newBranch
                if newNode != node: # Make sure to not to include our original node

                    newTreeList.append(getRootedSubtree(tree, newNode, newBranch)) # Recursively calls on the right and left subtree

    # Combine node name and branch length of tree with right and left subtrees
    newTree = ((nName, bLength),
                 newTreeList[0],
                 newTreeList[1])

    return newTree





def optimizeOneBranchHelper(leftSubtree: tuple, rightSubtree: tuple, seqD: dict, delta: float, branchlength : float) -> tuple[float, float]: # returns optimal branch length
    """ Returns optimal branch length given right and left subtree. Increments by delta every time it is run. 
        Helper for optimizeOneBranch
    INPUTS:
        • leftSubtree =  tuple, a rooted subtree in newick format on the "left" side of our branch
        • rightSubtree = tuple, a rooted subtree in newick format on the "right" side of our branch
        • seqD = dict, a dictionary of the node sequences
        • delta = float, an increment to multiply or divide by the current branch length, eg 1.01
        • branchlength = the length of the current branch
    OUTPUTS:
        • a float, the local optimal branch length of the current branch
        • a float, the likelihood of the tree with the local optimal branch length    
    """

    # Get score of current tree from with current branch length
    currentScore = branchLenLikelihood(leftSubtree, rightSubtree, seqD, branchlength)
    currentTuple = (currentScore, branchlength)
    
    # Get score of tree with increased branch length by delta
    timesDelBranchLength = branchlength * delta # multiply by delta
    timesDelScore = branchLenLikelihood(leftSubtree, rightSubtree, seqD, timesDelBranchLength)
    timesDelTuple = (timesDelScore, timesDelBranchLength)

    # Get score of the tree with decreased branch length by delta
    divDelBranchLength = branchlength / delta # divide by delta
    divDelScore = branchLenLikelihood(leftSubtree, rightSubtree, seqD, divDelBranchLength)
    divDelTuple = (divDelScore, divDelBranchLength)

    # Find which candidate length is the best
    maxTuple = max(currentTuple, timesDelTuple, divDelTuple)
    maxScore = maxTuple[0] # get max likelihood score
    bestBranchLength = maxTuple[1] # get best branch length
    
    if maxScore == currentScore: # Base case, at optimum where tree is not improving
        return [bestBranchLength, maxScore]
    else:
        assert(maxScore > currentScore) # If this is not true, something is breaking
        return optimizeOneBranchHelper(leftSubtree, rightSubtree, seqD, delta, bestBranchLength) # Recurse






def optimizeOneBranch(tree:Utree, seqD: dict, delta: float, branch :int) -> tuple[Utree,  float]:
    """ Returns optimal branch length and tree given right and left subtree. Increments by delta every time it is run.
    INPUTS:
        • tree =  Utree object, an unrooted tree
        • seqD = dict, a dictionary of the node sequences
        • delta = float, an increment to multiply or divide by the current branch length when optimizing, eg 1.01
        • branch = int, the index of branch to optimize
    OUTPUTS:
        • a Utree object, the unrotted tree with our optimized branch length
        • a float, the likelihood of the tree with the local optimal branch length    
    """
    leftNode, rightNode = tree.getNodes(branch) # gets nodes next to inputted branch
    leftSubtree = getRootedSubtree(tree, leftNode, branch) # get rooted version of left subtree
    rightSubtree = getRootedSubtree(tree, rightNode, branch) # get rooted version of right subtree
    currentBranchLength = tree.getBranchLength(branch) # retrieves branch length from Utree object

    # Use helper function to recursively optimize branch length and get likelihood score
    newBestBranchLength, newBestScore = optimizeOneBranchHelper(leftSubtree, rightSubtree, seqD, delta, currentBranchLength)

    tree.setBranchLength(branch, newBestBranchLength) # Updates Utree object

    return (tree, newBestScore)






def optimizeBranches(tree:Utree, seqD:dict, delta:float) -> tuple[Utree,float]:
    """ Optimizes Utree objects with sequence in seqD optimizes with increments of delta
    INPUTS:
        • tree =  Utree object, an unrooted tree
        • seqD = dict, a dictionary of the sequences in tree
        • delta = float, an increment to multiply or divide by the current branch length when optimizing, eg 1.01
        • branch = int, the index of branch to optimize
    OUTPUTS:
        • a Utree object, the updated unrooted tree
        • a float, the likelihood of the tree with optimized lengths    
    """

    branches = tree.getAllBranches() # get branch indices to iterate through
    oldTreeScores = [0 for x in branches] # set array to store likelihood after each branch optimization of previous cycle
    improvedBools = [True for x in branches] # set array to record if likelihood improved from cycle to cycle

    while any(improvedBools): # Loop terminates when no tree has improved from the last cycle to this one
        treeScores = []

        # iterates through each branch optimizing the tree, saves score to treeScores
        for thisBranch in branches:
            tree, newScore = optimizeOneBranch(tree, seqD, delta, thisBranch)
            treeScores.append(newScore) 

        # Records the trees that have improved from last cycle with their corresponding branch, used to terminate while loop
        improvedBools = [treeScores[i] > oldTreeScores[i] for i in range(len(branches))]
        
        oldTreeScores = treeScores # updates oldTreeScores to current scores
    
    return (tree, newScore)