# bio185d, Harvey Mudd College
# Eliot Bush, Aug 2016
# Joe Wirth, Aug 2022

import math
from UnrootedTree import Utree

##### GLOBAL CONSTANTS #####
CHARS = ["A", "T", "C", "G"]


##### FILE IO #####
def loadFasta(fastaFN:str) -> dict:
    """ loadFasta:
            Accepts the filename of a fasta file as input. Reads the sequence
            data into a dictionary whose keys are the taxon names and whose va-
            lues are the corresponding sequences as strings. Returns the dict.
    """
    # initialize the output dictionary
    outD = dict()

    # open the fasta file
    fh = open(fastaFN, 'r')

    # for each line in the file
    for line in fh:
        # if the first character is a '>'
        if line[0] == ">":
            # extract the name (drop the '>' and '\n' from the string)
            seqName = line[1:-1]

            # set the value for the name to an empty string
            outD[seqName] = ""
        
        # otherwise, the line is sequence data
        else:
            # chop off a new line character if present
            if line[-1] == "\n":
                line = line[:-1]

            # append the sequence (all caps) to the growing string
            outD[seqName] += line.upper()

    # close the file
    fh.close()

    return outD
 

def loadTree(treeFN:str) -> Utree:
    """ loadTree:
            Accepts the filename for a tree in newick format. Constructs and
            returns the corresponding Utree object.
    """
    fh = open(treeFN, 'r')
    tree = Utree(fh.read())
    fh.close()
    return tree


def writeUtree(tree:Utree, outFN:str) -> None:
    """ writeUtree:
            Accepts a Utree object and a filename for saving. Writes the Utree
            to file in newick format. Does not return.
    """
    # open file, write newick string, close file
    fh = open(outFN, 'w')
    fh.write(tree.toNewick() + "\n")
    fh.close()


##### PRELIMINARY BRANCH LENGTHS #####
def preliminaryBranchLengths(tree:tuple, seqD:dict) -> tuple:
    """ preliminaryBranchLengths:
            Accepts a tree in tuple format ((root,dist),left,right) and a dict-
            ionary whose keys are tips in the tree and whose values are aligned
            sequences. Returns a tree in tuple format with the same topology as
            the input tree but with some preliminary branch lengths provided.
    """
    # determine the order the nodes should be merged to create the input tree
    mergeOrderL = __getMergeOrder(tree)

    # convert the alignment to a distance matrix
    distD = __makeDistanceMatrix(seqD)

    # get a list of all the tip names
    namesL = list(seqD.keys())
    
    # initialize a dictionary to keep track of the branch lengths of each node
    brLenD = dict()

    # for each pair of nodes to be merged (except the final node)
    for nameA,nameB in mergeOrderL[:-1]:
        # extract the subtrees for the specified pair from the input tree
        subA = __findMatchingSubTree(nameA, tree)
        subB = __findMatchingSubTree(nameB, tree) 

        # get the branch lengths for the the two nodes
        lenA,lenB=__getBranchLength(nameA, nameB, namesL, distD)

        # save the branch lengths in the dictionary
        brLenD[nameA] = lenA
        brLenD[nameB] = lenB

        # merge the two subtrees together
        mergedNode = __mergeNodes(subA, subB, lenA, lenB, tree)

        # extract the name of the newly merged node
        mergedName = mergedNode[0][0]
        
        # add the distances to the distance matrix
        __updateDistances(nameA, nameB, mergedName, namesL, distD)

        # remove the processed names from the list
        namesL.remove(nameA)
        namesL.remove(nameB)

        # add the newly merged node name to the list
        namesL.append(mergedName)

    # extract the names of the last two nodes to be merged
    nameA,nameB=namesL

    # the distance for this final branch will be split evenly
    finalDist = distD[nameA, nameB] / 2

    # add the final distance to the branch length dictionary
    brLenD[nameA] = finalDist
    brLenD[nameB] = finalDist
    
    # use the dictionary and the original tree to construct the final tree
    return __buildFinalTree(tree,brLenD)


##### LIKELIHOOD FUNCTIONS #####
def jukesCantorProbability(branchLength:float) -> float:
    """ jukesCantorProbability:
            Accepts a branch length as input. Returns the probability that a
            single site will differ according to the Jukes-Cantor substitution
            model.
    """
    # convert a length of 'None' to 0
    if branchLength is None:
        branchLength = 0

    return (3/4) * (1 - math.exp(-(4/3) * branchLength))


def scoreOneChar(tree:tuple, character:str, pos:int, seqD:dict, \
                                                           memo:dict) -> float:
    """ scoreOneCharacter:
            Accepts a tree in tuple format ((root,dist),left,right), a single
            character (DNA), a position indicating the index of the alignment,
            a dictionary containing sequence data, and a memo dictionary as in-
            puts. Calculates and returns the likelihood for that character to
            occur given the provided tree.
    """
    # parse the tree into its individual components
    (root,dist),left,right = tree
    
    # check the memo for the solution before proceeding
    if (tree, character, pos) in memo:
        return memo[tree, character, pos]
    
    # base case: input tree is a tip
    elif len(left) == 0:
        # if the characters match, then it is 100% likely to have occurred
        if seqD[root][pos] == character:
            likelihood = 1.0
        
        # if the character's don't match, then it can occur with 0% likelihood
        else:
            likelihood = 0.0
    
    # recursive case
    else:
        # initialize the likelihood for the data
        likelihood = 0.0

        # extract the lengths of the branches leading to the left/right subtree
        leftBranchLen = left[0][1]
        rightBranchLen = right[0][1]

        # calculate the probability of base changes along left/right subtrees
        probLeftChanged = jukesCantorProbability(leftBranchLen)
        probRightChanged = jukesCantorProbability(rightBranchLen)

        # for each base pair on the left subtree
        for leftChar in CHARS:
            # recurse along the left subtree with the character
            leftLikelihood = scoreOneChar(left,leftChar,pos,seqD,memo)

            # if the input character matches the left character
            if character == leftChar:
                # then multiply by the likelihood that no change occurred
                leftLikelihood *= 1 - probLeftChanged

            # if the input character does not match the left character
            else:
                # then multiply by the likelihood that it changed to leftChar
                leftLikelihood *= probLeftChanged / 3
            
            # for each base pair on the right subtree
            for rightChar in CHARS:
                # recurse along the right subtree with the character
                rightLikelihood = scoreOneChar(right,rightChar,pos,seqD,memo)
                
                # if the input character matches the right character
                if character == rightChar:
                    # then multiply by the likelihood that no change occurred
                    rightLikelihood *= 1 - probRightChanged
                
                # if the input character does not match the right character
                else:
                    # then multiply by likelihood that it changed to rightChar
                    rightLikelihood *= probRightChanged / 3
                
                # add the score for this combination of characters to the sum
                likelihood += leftLikelihood * rightLikelihood
    
    # save the result for this calculation in the memo
    memo[tree, character, pos] = likelihood
    
    return likelihood


def scoreOnePosition(tree:tuple, pos:int, seqD:dict) -> float:
    """ scoreOnePosition:
            Accepts a tree in tuple format ((root,dist),left,right), a positi-
            on indicating the index of the alignment, and a dictionary contain-
            ing sequence data as inputs. Calculates the likelihood for the data
            at a single position and returns it.
    """
    # initialize the likelihood
    likely = 0.0

    # for each possible base pair
    for char in CHARS:
        # add its likelihood to the sum weighted by the number of characters
        likely += scoreOneChar(tree, char, pos, seqD, {}) / len(CHARS)
    
    # average the likelihood across the number of bases
    return likely


def calculateLikelihood(tree:tuple, seqD:dict) -> float:
    """ calculateLikelihood:
            Accepts a tree in tuple format ((root,dist),left,right) and a dict-
            ionary whose keys are the names of the tips and whose values are
            the aligned sequences. Returns the likelihood of the tree given the
            alignment.
    """
    # determine the length of the alignment
    seqLen = len(list(seqD.values())[0])

    # initialize the likelihood 
    likelihood = 1

    # for each position in the alignment
    for position in range(seqLen):
        # calculate its likelihood and multiply it by the growing total
        likelihood *= scoreOnePosition(tree, position, seqD)
    
    return likelihood



##### PRIVATE HELPER FUNCTIONS #####
def __propDiff(seqA:str, seqB:str) -> float:
    """ propDiff:
            Accepts two aligned sequences (str) as inputs. Calculates the prop-
            ortion of differences (ignoring gaps). Returns the value.
    """
    # constant
    GAP = "-"

    # initialize variables for counting
    seqLenNoGap = 0
    numMismatch = 0
    
    # for each position in the alignment
    for idx in range(len(seqA)):
        # if neither of the characters are a gap
        if seqA[idx] != GAP and seqB[idx] != GAP:
            # increment the sequence length without gaps counter
            seqLenNoGap += 1

            # if the characters don't match
            if seqA[idx] != seqB[idx]:
                # increment the mismatch counter
                numMismatch += 1
    
    return numMismatch / seqLenNoGap


def __jukesCantor(propDiff:float) -> float:
    """ jukesCantor:
            Accepts a a proportion of observed differences as input. Applies
            the Jukes-Cantor method to calculate the "actual" number of changes
            and returns the value.
    """
    # input value should not be negative
    if propDiff < 0:
        raise ValueError("negative values are invalid")
    
    # input values cannot exceed 0.75
    if propDiff >= 0.75:
        raise ValueError("values ≥ 0.75 are invalid")
    
    return -(3/4) * math.log(1 - (4/3) * propDiff)


def __makeDistanceMatrix(seqD:dict) -> dict:
    """ makeDistanceMatrix:
            Accepts a dictionary whose keys are the names of the tips and whose
            values are the aligned sequences. Constructs a dictionary whose
            keys are pairs of names and whose values are the distances between
            them based on the Jukes-Cantor corrected distances of the aligned
            sequences. Returns the dictionary.
    """
    # initialize the output dictionary
    distD = dict()

    # for each pairwise combination of aligned sequences
    for sp1 in seqD.keys():
        for sp2 in seqD.keys():
            # calculate the jukes-cantor corrected distance bewteen them
            dist = __jukesCantor(__propDiff(seqD[sp1], seqD[sp2]))

            # store the values in the dictionary (forward/reverse comparisons)
            distD[sp1,sp2] = dist
            distD[sp2,sp1] = dist
    
    return distD


def __getMergeOrder(tree:tuple) -> list:
    """ getMergeOrder:
            Accepts a tree in tuple format ((root,dist),left,right) as input.
            Constructs and returns a list which shows the order in which the
            nodes should be merged in order to produce the input tree, starting
            from the tips.
    """
    # parse the tree into its individual components
    (root,dist),left,right = tree

    # initialize the output list
    outL = []

    # if tree is not a tip
    if left != ():
        # recurse on the left and right subtrees
        leftL = __getMergeOrder(left)
        rightL = __getMergeOrder(right)

        # save the results of the recursive calls
        outL.extend(leftL)
        outL.extend(rightL)

        # save the names of the left and right subtrees for the current node
        leftName = left[0][0]
        rightName = right[0][0]
        outL.append((leftName,rightName))

    return outL


def __findMatchingSubTree(name:str, tree:tuple) -> tuple:
    """ findMatchingSubTree:
            Accepts the name of an internal or external node (str) and a tree
            in tuple format ((root,dist),left,right) as inputs. Recursively na-
            vigates the tree until the specified internal node is found. Retur-
            ns the subtree as a tuple or False if a matching subtree was not
            found.
    """
    # parse the tree into its individual components
    (root,dist),left,right = tree

    # return the input if it is the subtree we are searching for
    if root == name:
        return tree
    
    # if we are at a tip, then a subtree was not found
    elif left == ():
        return False
    
    # otherwise, recurse on the left and right subtrees
    else:
        leftSearch = __findMatchingSubTree(name, left)
        rightSearch = __findMatchingSubTree(name, right)

        # if the left tree contains the desired node, then return it.
        if leftSearch:
            return leftSearch
        
        # otherwise, return the result from the right subtree
        return rightSearch


def __nodeSep(name:str, nodeL:list, distD:dict) -> float:
    """ nodeSep:
            Accepts the name of a node (str), a list of all unprocessed nodes,
            and a distance matrix (dict) as inputs. Calculates a measure of se-
            paration between the specified nodes and all other nodes. Returns
            the value as a float.
    """
    # initialize the output value
    distSum = 0

    # for each node in the node list
    for otherName in nodeL:
        # add the distance to the running sum; avoid self-vs-self comparisons
        if name != otherName:
            distSum += distD[(name,otherName)]

    # average the sum by number of nodes minus 2
    return(float(distSum)/(len(nodeL)-2))


def __getBranchLength(nameA:str, nameB:str, nodeL:list, distD:dict) -> \
                                                            tuple[float,float]:
    """ getBranchLength:
            Accepts two node names (str) that are connected, a list of unproce-
            ssed nodes, and a distance matrix (dict) as inputs. Calculates the
            branch lengths from the parental node and returns them as a tuple.

    """
    # extract the distance between the two input nodes
    dist = distD[(nameA,nameB)]

    # get a measure of the separation from each of the other nodes
    sepA = __nodeSep(nameA,nodeL,distD)
    sepB = __nodeSep(nameB,nodeL,distD)

    # calculate the branch lengths
    branchA = (dist + (sepA - sepB)) / 2
    branchB = (dist + (sepB - sepA)) / 2

    return(branchA,branchB)


def __findParentName(nodeA:tuple, nodeB:tuple, tree:tuple) -> str:
    """ findParentName:
            Accepts three trees in tuple format ((root,dist),left,right) where
            the first two trees are subtrees within the third tree. Recursively
            navigates the first tree until the parental node of the first two
            trees is found. Returns the name of the parental node or False if a
            parental node could not be found.
    """
    # parse the tree into its individual components
    (root,dist),left,right = tree

    # if the tree is a tip, then it is not possible to search further
    if left == ():
        return False

    # if the left and right subtrees match the two specified nodes
    elif (nodeA,nodeB) == (left,right) or (nodeA,nodeB) == (right,left):
        # then this is the parent we are looking for, return its name
        return root
    
    else:
        # recurse along the left and right subtrees
        leftSearch = __findParentName(nodeA,nodeB,left)
        rightSearch = __findParentName(nodeA,nodeB,right)

        # return the answer from the left subtree if one was found
        if leftSearch:
            return leftSearch
    
        # otherwise return the answer from the right subtree
        return rightSearch


def __mergeNodes(nodeA:tuple, nodeB:tuple, branchLenA:float, \
                                        branchLenB:float, tree:tuple) -> tuple:
    """ mergeNodes:
            Accepts two subtrees in tuple format ((root,dist),left,right) to be
            merged, the lengths of the branches leading to the subtrees, and
            the full tree in tuple format as inputs. Constructs a new node from
            the inputs and returns it.
    """
    # parse the subtrees into their individual components
    ((rootA,distA),leftA,rightA) = nodeA
    ((rootB,distB),leftB,rightB) = nodeB

    # construct the new subtrees with the adjusted branch lengths
    newNodeA=((rootA,branchLenA),leftA,rightA)
    newNodeB=((rootB,branchLenB),leftB,rightB)

    # determine what the root of this merged node should be named
    parentName = __findParentName(nodeA,nodeB,tree)

    # construct and return the merged node
    return ((parentName,0),newNodeA,newNodeB)


def __updateDistances(nameA:str, nameB:str, parentName:str, nodeL:list, \
                                                           distD:dict) -> None:
    """ updateDistances:
            Accepts two names of sibling nodes, the name of the siblings' pare-
            ntal node, a list of the unprocessed nodes, and a distance matrix
            (dict) as inputs. Calculates the distances between the parent node
            and all other nodes in nodeL and stores them in the provided dict.
            Does not return.
    """
    # for each node
    for curNode in nodeL:
        # if the node is not one of the sibling nodes
        if curNode != nameA and curNode != nameB:
            # calculate the total distance between both siblings and the node
            newDist = distD[nameA,curNode] + distD[nameB,curNode]

            # subtract the distance between the two sibling nodes
            newDist -= distD[nameA,nameB]

            # halve the distance; current value is twice the actual distance
            newDist /= 2

            # store the new distance in the dictionary
            distD[(parentName,curNode)] = newDist
            distD[(curNode,parentName)] = newDist


def __buildFinalTree(tree:tuple, brLenD:dict) -> tuple:
    """ buildFinalTree:
            Accepts a tree in tuple format ((root,dist),left,right) and a dict-
            ionary whose keys are the names of nodes in the tree and whose vals
            are the corresponding branch lengths for the nodes. Constructs and
            returns a new tree in tuple format with the distances specified in
            the input dictionary.
    """
    # parse the tree into its individual components
    (root,dist),left,right = tree

    # if the name of the tree is in the dictionary
    if root in brLenD.keys():
        # then replace its distance with the corresponding value
        dist = brLenD[root]

    # otherwise, set it to zero (this handles 'anc0' being absent)
    else:
        dist = 0
    
    # if the tree is not a tip, then recurse along the left/right subtrees
    if left != ():
        left = __buildFinalTree(left,brLenD)
        right = __buildFinalTree(right,brLenD)

    # construct and return the final tree
    return ((root,dist), left, right)

