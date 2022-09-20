from branches import *


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
