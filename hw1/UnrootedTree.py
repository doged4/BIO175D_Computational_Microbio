# UnrootedTree class
# bio 185d
# Eliot Bush, Aug 2016
# Joe Wirth, Aug 2022

from __future__ import annotations
import copy, re

class Utree:
    ##### CONSTRUCTOR #####
    def __init__(self,nwkStr:str) -> Utree:
        '''Creates a Utree object from a newick string'''

        # initialize member variables
        self.__nodeNameD:dict = dict()
        self.__nodeBranchD:dict = dict()
        self.__branchNodeD:dict = dict()
        self.__branchLenD:dict = dict()

        # read the newick string into the member variables
        self.__readNewickStr(nwkStr)

        # nodeBranchD and nodeNameD should have the same keys
        if sorted(self.__nodeBranchD.keys()) != sorted(self.__nodeNameD.keys()):
            raise ValueError("nodeBranchD and nodeNameD don't have the same keys.")

        # tip nodes should have one branch, internal nodes three
        for key in self.__nodeBranchD.keys():
            # find internal nodes
            if not self.isLeaf(key):
                if len(self.__nodeBranchD[key]) != 3:
                    raise ValueError("Internal nodes in nodeBranchD should be attached to 3 branches.")
            else:
                if len(self.__nodeBranchD[key]) != 1:
                    raise ValueError("Tips in nodeBranchD should be attached to 1 branch.")

        # every branch should have two nodes
        if not all(len(self.__branchNodeD[branch])==2 for branch in self.__branchNodeD):
            raise ValueError("Error in branchNodeD, each branch should be attached to two nodes.")

        # save the collection of branches
        branchList = list(self.__branchLenD.keys())
        branchList.sort()
        self.__branches = tuple(branchList)

        # save the collection of nodes
        nodeList = list(self.__nodeBranchD.keys())
        nodeList.sort()
        self.__nodes = tuple(nodeList)


    ##### OVERLOADS #####
    def __copy__(self, memo) -> Utree:
        """ allows copy.deepcopy(Utree) to work
        """
        # make a new placeholder tree
        out = Utree("(A,B);")

        # copy all member variables into the new object
        out.__nodeNameD = copy.deepcopy(self.__nodeNameD)
        out.__nodeBranchD = copy.deepcopy(self.__nodeBranchD)
        out.__branchNodeD = copy.deepcopy(self.__branchNodeD)
        out.__branchLenD = copy.deepcopy(self.__branchLenD)
        out.__branches = copy.deepcopy(self.__branches)
        out.__nodes = copy.deepcopy(self.__nodes)
        
        return out


    def __repr__(self) -> str:
        """ allows repr(Utree) to work
        """
        return str(self)
        

    def __str__(self) -> str:
        """ allows str(Utree) to work
        """
        # initialize a list and add the first two rows to it
        rowsL=[]
        rowsL.append(['Sp.','Nnum)', 'branch','Length','(Nnum','Sp.'])
        rowsL.append(['---','-----', '------','------','-----','---'])

        # for each branch
        for branch in self.__branchLenD.keys():
            # extract the two nodes that it is connected to
            node1,node2 = self.__branchNodeD[branch]

            # get the string for the first node
            nm1Str = ""
            if not re.match(r"^anc\d+$", self.getName(node1)):
                nm1Str = self.getName(node1)
            
            # get the string for the second node
            nm2Str = ""
            if not re.match(r"^anc\d+$", self.getName(node2)):
                nm2Str = self.getName(node2)

            # get the string for the branch length
            if self.__branchLenD[branch] == None:
                brLenStr = str(None)
            else:
                brLenStr = format(self.__branchLenD[branch],".4f")
            
            # add the current branch's data to the row
            rowsL.append([nm1Str,
                          str(node1) + ")",
                          str(branch),
                          brLenStr,
                          "(" + str(node2),
                          nm2Str])
        
        return Utree.__tableString(rowsL)


    def __eq__(self, other:Utree) -> bool:
        # TODO: implement this function

        # can only compare Utrees to other Utrees
        if type(other) is not Utree:
            raise ValueError("cannot compare Utree to " + str(type(other)))
        
        # fail by default. REMOVE THIS LINE BEFORE IMPLEMENTING!
        raise RuntimeError("an equality function has not been implemented")
    

    def __ne__(self, other:Utree) -> bool:
        """ returns a boolean indicating if two Utree objects do not represent
            equivalent trees (eg. (A,(B,C)) != (D,(A,B))).
        """
        return not self == other


    ##### PUBLIC METHODS #####
    def getAllBranches(self) -> tuple:
        """ getAllBranches:
                Accepts no inputs. Returns all of the branches as a tuple of
                integers.
        """
        return self.__branches


    def getAllNodes(self) -> tuple:
        """ getAllNodes:
                Accepts no inputs. Returns all of the nodes as a tuple of inte-
                gers.
        """
        return self.__nodes


    def getNodes(self, branch:int) -> tuple[int,int]:
        """ getNodes:
                Accepts a branch (int) as input. Returns the two nodes that the
                branch is connected to as a tuple.
        """
        return self.__branchNodeD[branch]


    def getBranches(self, node:int) -> tuple:
        """ getBranches:
                Accepts a node (int) as input. Returns the branches (one for 
                leaves; three for internal nodes) as a tuple of integers.
        """
        return self.__nodeBranchD[node]


    def getName(self, node:int) -> str:
        """ getName:
                Accepts a node (int) as input. Returns the name for that node.
        """
        return self.__nodeNameD[node]


    def isLeaf(self, node:int) -> bool:
        """ isLeaf:
                Accepts a node number as input. Returns a boolean indicating
                whether or not the node is a leaf in the calling object.
        """
        # constant
        GREP_FIND = r"^anc\d+$" # finds internal nodes

        # finds and returns a boolean
        return not bool(re.match(GREP_FIND, self.getName(node)))


    def getBranchLength(self, branch:int) -> float:
        """ getBranchLength:
                Accepts a branch (int) as input. Returns its length as a float.
        """
        return self.__branchLenD[branch]


    def setBranchLength(self, branch:int, branchLen:float) -> None:
        """ setBranchLength:
                Accepts a branch (int) and a length. Updates the length of the
                specified branch within the calling object. Does not return.
        """
        if branchLen is not None:
            branchLen = float(branchLen)
        
        self.__branchLenD[branch] = branchLen


    def printNodeNameSeq(self,dataD:dict) -> None:
        '''Given a dict of sequences, print node number, name, seq for each node.'''
        print("Node number, name and sequence for tips:")
        rowL=[]
        for node in self.__nodes:
            if self.isLeaf(node):
                rowL.append([str(node),self.getName(node),dataD[self.getName(node)]])

        print(Utree.__tableString(rowL))


    def toNewick(self) -> str:
        """ toNewick:
                Accepts no inputs. Generates and returns a string in newick
                format corresponding to the calling object. The resulting tree
                is arbitrarily rooted.
        """
        treeT = self.toTuple()
        return Utree.__toNewickHelper(treeT) + ";"


    def toTuple(self) -> tuple:
        """toTuple:
                Accepts no inputs. Uses a recursive helper function to generate
                the tree in tuple format ((root,dist),left,right) that corresp-
                onds to the calling object. Returns the tuple. The resulting
                tree is arbitrarily rooted.
        """
        # constants
        ROOT_BRANCH = 1
        ROOT_NAME = "anc0"
        DEFAULT_ROOT_DIST = float(0)
        
        # this is an unrooted tree; must process the "root" differently
        # get the nodes that are connected to the root branch
        leftNode,rightNode = self.getNodes(ROOT_BRANCH)

        # if there is a branch length for the root branch
        if self.getBranchLength(ROOT_BRANCH) is not None:
            # split the distance of the root branch between the two subtrees
            dist = self.getBranchLength(ROOT_BRANCH) / 2
        
        # otherwise the distance is None
        else:
            dist = None

        # use recursive helper function to traverse left and right subtrees
        leftSub = self.__toTupleHelper(leftNode, dist, {ROOT_BRANCH})
        rightSub = self.__toTupleHelper(rightNode, dist, {ROOT_BRANCH})

        # if the branch length is None
        if dist is None:
            # then the root's branch length should also be None
            rootDist = None
        
        # otherwise the root's branch length should be the default value
        else:
            rootDist = DEFAULT_ROOT_DIST

        # construct and return the tuple
        return ((ROOT_NAME,rootDist), leftSub, rightSub)


    def fromTuple(tree:tuple) -> Utree:
        """ fromTuple:
                Accepts a tree in tuple format. Constructs a newick string from
                it and uses that to build and return a Utree object.
        """
        nwkStr = Utree.__toNewickHelper(tree)
        return Utree(nwkStr)


    ##### PRIVATE METHODS #####
    def __readNewickStr(self, nwkStr:str) -> None:
        """ readNewickStr:
                Accepts a newick string as input. Extracts the relevant data
                from the string and stores the data in the calling object's 
                member variables. Does not return.
        """
        # constants
        FIND_DIST = r":(\d[^,&^\)]{0,})" # finds the distance
        FIND_TAXA = r"([\(|,]+)([^,&^\)]+)" # finds taxa
        FIND_SPACE = r"([,|\(|\)]) " # finds whitespace after comma or paren

        # chop off newlines and semicolons
        while nwkStr[-1] in ["\n", ";"]:
            nwkStr = nwkStr[:-1]
        
        # initialize variables
        branchLenD = dict()
        
        Utree.__getAncestorDistances(nwkStr, branchLenD)

        ### make a tuple of only the taxa
        # remove the distances from the string
        noDistStr = re.sub(FIND_DIST, r"", nwkStr)

        # remove any trailing spaces
        noDistStr = re.sub(FIND_SPACE, r"\1", noDistStr)
        
        # add quotes around the taxa
        noDistStr = re.sub(FIND_TAXA, r'\1"\2"', noDistStr)

        # convert the string to a tuple
        noDistT = eval(noDistStr)

        # call recursive helper function to get tip distances
        Utree.__getTipDistances(noDistT, nwkStr, branchLenD)

        # call recursive helper function to build the tuple tree
        treeT = Utree.__buildTupleTree(noDistT, branchLenD, {"anc":0})

        # populate member variables
        self.__getDataFromTupleTree(treeT)


    def __getParenthesesIndices(nwkStr:str) -> dict:
        """ getParenthesesIndices:
                Accepts a newick string as input. Uses a stack to identify the
                indices for each pair of open/close parentheses. Returns a dict
                whose keys are the names of the internal nodes and whose values
                are the indices (tuple) required to extract the node's string
                from the input string.
        """
        # initialize variables for looping
        idxOpenL = list()  # the stack used to track opened paren
        ancL = list()  # the stack used to track the name of the internal nodes
        out = dict()  # the dictionary to be returned
        curAnc = 0

        # for each index and corresponding character in the string ...
        for idx,char in enumerate(nwkStr):
            # if the character is an open paren
            if char == "(":
                # add the index and the next unused ancestor to the stacks
                idxOpenL.append(idx)
                ancL.append("anc" + str(curAnc))

                curAnc += 1
            
            # if the character is a closed parent
            elif char == ")":
                # get the most recently added ancestor and index in the stacks
                try:
                    out[ancL.pop()] = (idxOpenL.pop(), idx + 1)
                
                # if nothing in the stacks, then there is an unmatched paren
                except:
                    raise ValueError("unmatched parentheses") 
        
        # if items remain in the stack after looping, then unmatched paren
        if len(idxOpenL) > 0:
            raise ValueError("unmatched parentheses")
        
        return out
        

    def __getAncestorDistances(nwkStr:str, branchLenD:dict) -> None:
        """ getAncestorDistances:
                Accepts a newick string and a dictionary of branch lengths as
                inputs. Extracts the distances for the internal nodes and saves
                them in the provided dictionary. Does not return.
        """
        # constants
        FIND_FLOAT = r"(^[\d|\.]+).+$" # finds a float

        # get a dict of open/close parentheses locations in the string
        matchedParen = Utree.__getParenthesesIndices(nwkStr)

        # sort the list to ensure that the inner-most ancestors come first
        ancL = sorted(list(matchedParen.keys()), reverse=True)

        # for each internal node
        for anc in ancL:
            # extract the matching parentheses' indices from the dict
            openIdx,closeIdx = matchedParen[anc]

            # remove all text (including colon) before the distance value
            subStr = nwkStr[closeIdx+1:]

            # use regexp to extract the distance from the substring
            dist = re.sub(FIND_FLOAT, r"\1", subStr)
            
            # if no distance was found, then set it to None
            try:
                branchLenD[anc] = float(dist)
            except:
                branchLenD[anc] = None


    def __getTipDistances(nwkT:tuple, nwkStr:str, branchLenD:dict) -> None:
        """ getTipDistances:
                Accecpts a tuple (literal eval of distance-less newick string)
                the original newick string, and a dictionary of branch lengths
                as inputs. Recursively navigates the tuple to extract the dist-
                ances for each tip from the newick string and saves them in the
                provided dictionary. Does not return.
        """ 
        # get left and right sub trees
        left, right = nwkT

        # if the left subtree is a string, then it is a tip
        if type(left) is str:
            # extract the distance for the corresponding tip
            leftDist = Utree.__extractTipDistanceFromNewick(left, nwkStr)

            # save the tip's branch length
            branchLenD[left] = leftDist
        
        # otherwise, the left subtree is an internal node
        else:
            # recurse on the left subtree
            Utree.__getTipDistances(left, nwkStr, branchLenD)

        # if the right subtree is a string, then it is a tip
        if type(right) is str:
            # extract the distance for the corresponding tip
            branchLenD[right] = Utree.__extractTipDistanceFromNewick(right, nwkStr)

        else:
            # recurse on the right subtree
            Utree.__getTipDistances(right, nwkStr, branchLenD)


    def __extractTipDistanceFromNewick(tipName:str, nwkStr:str) -> float:
        """ extractTipDistanceFromNewick:
                Accepts two strings, the name of a tip and a newick string, as
                inputs. Uses regular expressions to extract the distance for
                the given tip's branch. Returns the distance as a float or as
                'None' if there is no distance present.
        """
        # constants
        GREP_PREFIX = r"^.+"
        GREP_SUFFIX = r":([^\(&^\)&^,]+)[\(|\)|,].{0,}$"
        GREP_REPL = r"\1"

        # construct a regex to find the distance for the tip
        find = GREP_PREFIX + tipName + GREP_SUFFIX

        # get the distance from the string
        dist = re.sub(find, GREP_REPL, nwkStr)

        # convert to float
        try:
            dist = float(dist)
        
        # set to None if it cannot be coerced to float
        except:
            dist = None
        
        return dist


    def __buildTupleTree(nwkT:tuple, branchLenD:dict, curAnc:dict) -> tuple:
        """ buildTupleTree:
                Accecpts a tuple (literal eval of distance-less newick string),
                a dictionary of the branch lengths, and a dictionary to track
                the most recently used ancestor number as inputs. Recursively
                builds a tuple tree '((root,dist),left,right)' and returns it.
        """
        # constant
        ANC = "anc"

        # parse data from the input tuple
        left,right = nwkT

        # get the current ancestor string and update the value in the dict
        anc = ANC + str(curAnc[ANC])
        curAnc[ANC] += 1

        # base case 1: left is a tip
        if type(left) is str:
            # look up the distance
            dist = branchLenD[left]

            # construct the left subtree
            leftSub = ((left,dist), (), ())
        
        # otherwise recurse on left subtree
        else:
            leftSub = Utree.__buildTupleTree(left, branchLenD, curAnc)

        # base case 2: right is tip
        if type(right) is str:
            # look up the distance
            dist = branchLenD[right]

            # construct the right subtree
            rightSub = ((right,dist), (), ())
        
        # otherwise recurse on right subtree
        else:
            rightSub = Utree.__buildTupleTree(right, branchLenD, curAnc)
        
        # compile and return the tuple
        return ((anc, branchLenD[anc]), leftSub, rightSub)
        
    
    def __getDataFromTupleTree(self, treeT:tuple) -> None:
        """ getDataFromTupleTree:
                Accepts a tree in tuple format ((root,dist),left,right) and
                uses it to populate the member variables of the calling object.
                Does not return.
        """
        # parse the tuple tree into its individual components
        (root,dist),left,right = treeT

        # initialize variables
        rootNode = 0
        tempRootBranch = 0 # an branch attached to this root node we'll remove later
        branchCount = 0
        nodeCount = 1

        # initialize a list for the root node
        self.__nodeBranchD[rootNode] = []
        
        # use the left subtree to populate the member variables
        branchCount,nodeCount = self.__populateMemberVariables(left, rootNode,
                                                        branchCount, nodeCount)

        # use the right subtree to populate the member variables
        branchCount,nodeCount = self.__populateMemberVariables(right, rootNode,
                                                        branchCount, nodeCount)


        # now edit these dictionaries to remove the root node and the
        # tempRootBranch branch attached to that root node. And update the
        # connections accordingly

        branch1,branch2 = self.__nodeBranchD[rootNode]
        branchToKeep = branch1 if branch2==tempRootBranch else branch2

        # must update branchToKeep entry in branchNodeD
        node1,node2 = self.__branchNodeD[tempRootBranch]
        nodeToReplaceWith = node1 if node2==rootNode else node2

        node1,node2 = self.__branchNodeD[branchToKeep]
        if node1 == rootNode:
            self.__branchNodeD[branchToKeep] = [nodeToReplaceWith,node2]
        else:
            self.__branchNodeD[branchToKeep] = [node1,nodeToReplaceWith]

        # if the branchToKeep has a length
        if self.__branchLenD[branchToKeep] is not None:
            # then add the length of tempRootBranch to branchToKeep's length
            self.__branchLenD[branchToKeep] +=self.__branchLenD[tempRootBranch]
        
        # update the nodeToReplaceWith entry in nodeBranchD. Take out
        # tempRootBranch and put in branchToKeep
        ind = self.__nodeBranchD[nodeToReplaceWith].index(tempRootBranch)
        self.__nodeBranchD[nodeToReplaceWith][ind] = branchToKeep

        # now get rid of rootNode and tempRootBranch entries in nodeBranchD,
        # branchNodeD, and branchLenD 
        del self.__nodeBranchD[rootNode]
        del self.__branchNodeD[tempRootBranch]
        del self.__branchLenD[tempRootBranch]

        # convert the values of branchNodeD from lists to tuples
        for key in self.__branchNodeD.keys():
            self.__branchNodeD[key] = tuple(self.__branchNodeD[key])
        
        # convert the values of nodeBranchD from lists to tuples
        for key in self.__nodeBranchD.keys():
            self.__nodeBranchD[key] = tuple(self.__nodeBranchD[key])


    def __populateMemberVariables(self, treeT:tuple, parentNode:int, \
                                      branchCount:int, nodeCount:int) -> tuple:
        """ populateMemberVariables:
                Accecpts a tree in tuple format ((root,dist),left,right), an
                integer indicating the parent node, an int indicating the curr-
                ent branch number, and an integer indicating the current node
                number. Recursively navigates the tuple to populate the calling
                object's member variables. Returns a tuple containing the upda-
                ted branch number and node number.
        """
        # store values from inputs and increment the counts
        thisBranch = branchCount
        thisNode = nodeCount
        branchCount += 1
        nodeCount += 1

        # parse tuple tree into its individual components
        (root,dist),left,right = treeT
        
        # populate the member variables with relevant data
        self.__branchNodeD[thisBranch] = [thisNode,parentNode]
        self.__branchLenD[thisBranch] = dist
        self.__nodeBranchD[thisNode] = [thisBranch]
        self.__nodeBranchD[parentNode].append(thisBranch)
        self.__nodeNameD[thisNode] = root

        # if the current node is not a tip
        if not left == ():
            # recurse along the left subtree
            branchCount,nodeCount = self.__populateMemberVariables(left,
                                                                   thisNode,
                                                                   branchCount,
                                                                   nodeCount)
            
            # recurse along the right subtree
            branchCount,nodeCount = self.__populateMemberVariables(right,
                                                                   thisNode,
                                                                   branchCount,
                                                                   nodeCount)

        return branchCount,nodeCount


    def __toTupleHelper(self,node:int,dist:float,branchesSeen:set) -> tuple:
        """ toTupleHelper:
                Accepts a node, a distance (branch len), and a set of branches
                that have already been seen as inputs. Recursively navigates
                the calling object to construct a tuple representation of the
                tree ((root,dist),left,right). Returns the tuple.
        """
        # get the name of the root node
        root = self.getName(node)

        # base case: the node is a tip
        if self.isLeaf(node):
            return ((root,dist), (), ())
        
        # recursive case
        else:
            # get a set of all the branches that the node is connected to
            connectedBranches = set(self.getBranches(node))

            # ignore the branches we have already seen
            connectedBranches.difference_update(branchesSeen)

            # add these new branches to the set
            branchesSeen.update(connectedBranches)

            # there will be exactly two unseen branches; left/right subtrees
            leftBranch,rightBranch = connectedBranches

            # get the distances for the left/right subtrees
            leftDist = self.getBranchLength(leftBranch)
            rightDist = self.getBranchLength(rightBranch)

            # get the node in the left subtree (ignore current node)
            leftNode = set(self.getNodes(leftBranch))
            leftNode.difference_update({node})
            leftNode = leftNode.pop()

            # get the node in the right subtree (ignore current node)
            rightNode = set(self.getNodes(rightBranch))
            rightNode.difference_update({node})
            rightNode = rightNode.pop()

            # recurse along the left subtree
            leftSub = self.__toTupleHelper(leftNode,leftDist,branchesSeen)

            # recurse along the right subtree
            rightSub = self.__toTupleHelper(rightNode,rightDist,branchesSeen)
        
            return ((root,dist),leftSub,rightSub)


    def __toNewickHelper(treeT:tuple) -> str:
        """ toNewickHelper:
                Accepts a tree in tuple format ((root,dist),left,right). Recur-
                sively navigates the tree to create a string corresponding to
                the tree in newick format. Returns the string.
        """
        # parse the tree into its individual components
        (root,dist),left,right = treeT

        # base case: at a tip
        if left == ():
            outStr = root
        
        # recursive case
        else:
            # '( + left + , + right + )'
            outStr = "("
            outStr += Utree.__toNewickHelper(left)
            outStr += ","
            outStr += Utree.__toNewickHelper(right)
            outStr += ")"

        # append the distance if it is specified
        if dist != 0 and dist is not None:
            outStr += ":" + str(dist)

        return outStr


    def __tableString(inputList:list, indent:int=2) -> str:
        """ tableString:
                Accepts a matrix of data stored as a list of lists whose sub-
                lists are rows and an integer indicating the number of spaces 
                to put in front of each row as inputs. Converts the list of
                lists to a nicely formatted string and returns it.
        """
        # get max width for each column
        outStr=''
        colMax=[]
        for col in range(len(inputList[0])):
            mx=0
            for row in inputList:
                if len(row[col]) > mx:
                    mx = len(row[col])
            colMax.append(mx)
        
        # print
        for row in inputList:
            for col in range(len(row)):
                row[col]=row[col]+' ' * (colMax[col]-len(row[col]))
                
        for row in inputList:
            outStr += " "*indent + " ".join(row).rstrip()+'\n'

        return outStr

          
    def __getAllTips(treeT:tuple) -> list:
        """ getAllTips:
                Accepts a tree in tuple format (root,left,right) as input. Ret-
                urns a list of all the tips in the tree.
        """
        # parse the tuple tree
        root,left,right = treeT

        # base case: at a tip; return the tip name as a list
        if left == ():
            return [root]
        
        # recursive case: internal node
        else:
            # get all the tips in the left sub tree
            outL = Utree.__getAllTips(left)

            # add all the tips from right sub tree
            outL.extend(Utree.__getAllTips(right))

            return outL


    def __treePathFinder(self, origin:int, end:int, previous:int=None) -> \
                                                                   list[tuple]:
        """ treePathFinder:
                Accepts an origin node (int), an end node (int) and (optional) 
                the previous node as inputs. Recursively navigates the tree to 
                create a list of paired nodes such that the list reveals the
                path between the origin and the end nodes. Returns the path as
                a list of tuples.
        """
        # get all the nodes connected to the origin
        allConnxns = self.__getAllConnections(origin)

        # base case: the end is connected to the origin
        if end in allConnxns:
            return [(origin,end)]
        
        # recursive case: iterate through the origin's connections
        outL = list()
        for connxn in allConnxns:
            # only process connections that aren't leaves or the previous node
            if not self.isLeaf(connxn) and connxn != previous:
                # recurse on connxn towards the specified end; specify previous
                path = self.__treePathFinder(connxn, end, previous=origin)
                
                # if a path was found ...
                if path != []:
                    # ... then append the origin and its connection to the list ...
                    outL.append((origin, connxn))
                    # ... and add the path that was found to the list
                    outL.extend(path)

        return outL


    def __getAllConnections(self, node:int) -> set[int]:
        """ getAllConnections:
                Accepts  a node (int) as input. Constructs a set of all the
                nodes (int) that are connected to the provided node. Returns
                the set.
        """
        # get all the branches connected to the origin
        branches = self.getBranches(node)

        # initialize a set of nodes connected to the origin
        allConnxns = set()

        # for each branch
        for branch in branches:
            # get the nodes connected to the branch
            nodes = self.getNodes(branch)

            # for each connected node
            for n in nodes:
                # add it to the set unless it is the origin itself
                if n != node:
                    allConnxns.add(n)
        
        return allConnxns


    def __sumPathDistance(self, path:list[tuple]) -> float:
        """ sumBranchLength:
                Accepts a list of paired nodes (tuple) representing the path 
                between two nodes as inputs. Calculates and returns the total 
                distance along the path of paired nodes as a float.
        """
        # initialize length
        totalLen = 0.0

        # for each node-pair in the list
        for tip1,tip2 in path:
            # identify the branch connecting the pair
            branches1 = set(self.getBranches(tip1))
            branches2 = set(self.getBranches(tip2))
            branch = branches1.intersection(branches2).pop()

            # get thr branch's length
            brLen = self.getBranchLength(branch)
            
            # handle None lengths
            if brLen is None:
                brLen = 0
            
            # update the total length
            totalLen += brLen
        
        return totalLen
    

    def __removeDistancesFromTuple(treeT:tuple) -> tuple:
        """ removeDistancesFromTuple:
                Accepts a tree in tuple format ((root,dist),left,right) and
                constructs and returns a tree in tuple format without dist-
                ances (root,left,right).
        """
        # parse tuple tree
        (root,dist),left,right = treeT

        # if the tree is not a leaf, then recurse on the left and right 
        if left != ():
            left = Utree.__removeDistancesFromTuple(left)
            right = Utree.__removeDistancesFromTuple(right)
        
        # return the tree without distances
        return root,left,right


    def __reroot(treeT:tuple, tip:str) -> tuple:
        """ reroot:
                Accepts a tree in tuple format (root,left,right) WITHOUT dist-
                ances and a tip name (str) as inputs. Returns a tree that is
                rooted on the specified tip.
        """
        # constants
        ERR_MSG = 'input tree must be in tuple format WITHOUT distances'

        # make sure input is a tuple
        if type(treeT) is not tuple:
            raise ValueError(ERR_MSG)
        
        # check for distances
        try:
            # fail if they are present
            (root,dist),left,right = treeT
            raise ValueError(ERR_MSG)
        except: pass

        # handle if the spcified tip is not in the tree
        if tip not in Utree.__getAllTips(treeT):
            raise ValueError("specified tip is not a leaf")

        # parse the tree
        root, left, right = treeT

        # continue to rotate the tree until the specified tip is the root
        while (left[0] != tip) and (right[0] != tip):
            # if the tip is in the left subtree, then rotate it left
            if tip in Utree.__getAllTips(left):
                # parse the left subtree
                lName, A, B = left

                # if the tip is in the "left-left" subtree
                if tip in Utree.__getAllTips(A):
                    # ((A,B),right) => (A,(B,right))
                    left = A
                    right = (lName, B, right)
                
                # if the tip is in the "left-right" subtree
                else:
                    # ((A,B),right) => (B,(A,right))
                    left = B
                    right = (lName, A, right)
            
            # if the tip is in the right subtree, then rotate it right
            else:
                # parse the right subtree
                rName, C, D = right

                # if the tip is in the "right-left" subtree
                if tip in Utree.__getAllTips(C):
                    # (left,(C,D)) => ((left,D),C)
                    left = (rName, left, D)
                    right = C
                
                # if the tip is in the "right-right" subtree
                else:
                    # (left,(C,D)) => ((left,C),D)
                    left = (rName, left, C)
                    right = D
        
        return (root, left, right)


    def __getDistanceBetweenNodes(self, node1:int, node2:int) -> float:
        """ getDistanceBetweenNodes:
                Accepts two nodes (int) as inputs. Calculates and returns the
                distance between the two nodes as a float.
        """
        # get the path between the two nodes of interest
        pathL = self.__treePathFinder(node1, node2)

        # sum the distances along the path and return it
        return self.__sumPathDistance(pathL)

 