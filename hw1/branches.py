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
       • a float, the likelihood given the branch length, calculated by ((probchanged or same) * likelihood_right * likelihoodleft * 1/4)
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
    # TODO: comment, rename?
    """ Returns optimal branch length given right and left subtree. Increments by delta
    """

    currentScore = branchLenLikelihood(leftSubtree, rightSubtree, seqD, branchlength)
    currentTuple = (currentScore, branchlength)
    # print(currentScore)
    
    timesDelBranchLength = branchlength * delta
    timesDelScore = branchLenLikelihood(leftSubtree, rightSubtree, seqD, timesDelBranchLength)
    timesDelTuple = (timesDelScore, timesDelBranchLength)

    
    divDelBranchLength = branchlength / delta
    divDelScore = branchLenLikelihood(leftSubtree, rightSubtree, seqD, divDelBranchLength)
    divDelTuple = (divDelScore, divDelBranchLength)

    maxTuple = max(currentTuple, timesDelTuple, divDelTuple)
    maxScore = maxTuple[0]
    bestBranchLength = maxTuple[1]
    
    if maxScore == currentScore: # Base case, at optimum
        return [bestBranchLength, maxScore]
    else:
        assert(maxScore > currentScore)
        return optimizeOneBranchHelper(leftSubtree, rightSubtree, seqD, delta, bestBranchLength)






def optimizeOneBranch(tree:Utree, seqD: dict, delta: float, branch :int) -> tuple[Utree,  float]: # used to return new best branch as well
    # TODO: comment this function
    leftNode, rightNode = tree.getNodes(branch) # gets nodes next to inputted branch
    leftSubtree = getRootedSubtree(tree, leftNode, branch)
    rightSubtree = getRootedSubtree(tree, rightNode, branch)
    currentBranchLength = tree.getBranchLength(branch)

    newBestBranchLength, newBestScore = optimizeOneBranchHelper(leftSubtree, rightSubtree, seqD, delta, currentBranchLength)

    tree.setBranchLength(branch, newBestBranchLength)
    # return (tree, newBestBranchLength, newBestScore)
    return (tree, newBestScore)
    # timesDelBranchLength = currentBranchLength + delta
    # timesDelScore = branchLenLikelihood(leftSubtree, rightSubtree, seqD, timesDelBranchLength)

    # divDelBranchLength = currentBranchLength - delta
    # divDelScore = branchLenLikelihood(leftSubtree, rightSubtree, seqD, divDelBranchLength)

    


# def scoreThisTree(tree:Utree, seqD: dict, startBranch : int) -> int:
#     # TODO: delete this function, or integrate it better into optimizeOneBranchLength

#     leftNode, rightNode = tree.getNodes(startBranch)
#     leftSubtree = getRootedSubtree(tree, leftNode, startBranch)
#     rightSubtree = getRootedSubtree(tree, rightNode, startBranch)
#     currentBranchLength = tree.getBranchLength(startBranch)

#     currentScore = branchLenLikelihood(leftSubtree, rightSubtree, seqD, currentBranchLength)

#     return currentScore







def optimizeBranches(tree:Utree, seqD:dict, delta:float) -> tuple[Utree,float]:
    # TODO: write this function, comment, 
    newTree =  tree
    # newTree = deepcopy(tree)

    branches = newTree.getAllBranches()

    # GETS WORSE IN LOOP THROUGH BRANCHES??
    # cScore = scoreThisTree(newTree, seqD, 1)
    # print(cScore)
    # treeScores.append(cScore)
    oldMinScore = 0
    while True:
        treeScores = []
        for thisBranch in branches:
            tempnewTree, newScore = optimizeOneBranch(newTree, seqD, delta, thisBranch)
            # aScore = scoreThisTree(newTree, seqD, 1)
            # print(str(aScore) + " and " + str(newScore) + " : " + str(newScore - aScore))
            # assert (aScore == newScore)
            newTree = tempnewTree
            treeScores.append(newScore)

        assert oldMinScore <= treeScores[-1]
        # if not treeScores[-1] > treeScores[0] * (1 + 10**(-10)):
        if oldMinScore == treeScores[-1] :
            return (newTree, newScore)
        oldMinScore = treeScores[0]

    # print(treeScores) # They are getting worse over time??
    # assert treeScores[-1] > cScore # get rid
    # assert treeScores[-1] > cScore * (1 + 10**(-10))# get ri
    # print(scoreThisTree(tree, seqD, 1) - scoreThisTree(newTree, seqD, 1))
    # if  treeScores[-1] > cScore * (1 + 10**(-10)):

    


# ##### RUNNER for testing #####
# def test():
#     print('test1')
#     __main('example1.nwk','example1.fna',1.01,'example1Out.nwk')
#     print()
#     print('test2')
#     __main('example2.nwk','example2.fna',1.01,'example2Out.nwk')


# # WRAPPER FUNCTION
# def __main(treeFN:str, fastaFN:str, delta:float, outTreeFN:str) -> None:
#     """ main:
#             Accepts a tree filename, a fasta filename, a value to adjust branch
#             lengths, and a output filename as inputs. Optimizes the branch len-
#             gths for the provided tree. Prints the statistics of the tree and
#             saves it to the specified file. Does not return.
#     """
#     # constant
#     EOL = "\n"
    
#     # read data into memory
#     utree = loadTree(treeFN)
#     seqD = loadFasta(fastaFN)

#     # get some preliminary branch lengths
#     treeT = preliminaryBranchLengths(utree.toTuple(),seqD)
#     utree = Utree.fromTuple(treeT)
    
#     # optimize the tree
#     utree,likely = optimizeBranches(utree,seqD,delta)

#     # print the tree data for the final tree
#     print("Final branch lengths:")
#     print(utree, EOL)

#     print("Final Likelihood:", format(likely,".5e"))

#     # write the final tree to file
#     writeUtree(utree,outTreeFN)
    
#     print("Tree (arbitrarily rooted) with branch lengths written to:",outTreeFN)
    

# seqD = loadFasta('example1.fna')
# loadedtree = loadTree('example1.nwk')
# temptree  = preliminaryBranchLengths(loadedtree.toTuple(),seqD)
# ttree = Utree.fromTuple(temptree)
