# Joe Wirth
# Bio185d, Harvey Mudd College
# August 2022

import re

def findMotif(seqD:dict, motifRE:str) -> dict[str:list[re.Match]]:
    """"For motifRE regex string, returns dictionary of matched sequences with an array of all of their match objects for the sequence."""
    matchD = dict()
    # For each sample name in dictionary, do find to get iterable of match objects
    for name_key in seqD.keys():
        match_iter = re.finditer(motifRE, seqD[name_key])
        # List comprehension match iterable to array of match objects
        matches = [x for x in match_iter]
        
        # Only report matches if sample has at least one match
        if len(matches) != 0:
            matchD[name_key] = matches
    return matchD



def loadFasta(fastaFN:str) -> dict[str:str]:
    """ readFasta:
            Accepts a filename and returns a dictionary whose keys are the name
            of the sequence and whose values are the sequences themselves. Ret-
            urns the dictionary.
    """
    # initialize output dictionary
    seqD = dict()

    # open the file
    fh = open(fastaFN, 'r', encoding='utf-8')

    # for each line
    for line in fh.readlines():
        # truncate any newline characters
        if line[-1] == '\n':
            line = line[:-1]

        # grab new sequences
        if line[0] == ">":
            key = line[1:]
            seqD[key] = ""
        
        # populate sequences
        else:
            seqD[key] += line
    
    fh.close()
    return seqD



"""
Using findMotif.
Please put your answers to the questions within this docstring.
1. Ferredoxins
    regular expression: r'C..C..C.{3}CP'
    matching protein names (as a set): {'SPO_RS00735', 'SPO_RS01835', 'SPO_RS02970', 'SPO_RS04205', 'SPO_RS07240', 'SPO_RS07920', 'SPO_RS08845', 'SPO_RS09040', 'SPO_RS09130', 'SPO_RS09160', 'SPO_RS14070', 'SPO_RS17625', 'SPO_RS17845', 'SPO_RS18025', 'SPO_RS18785', 'SPO_RS19100']}
2. IQ-calmodulin-binding proteins
    regular expression: r'[FILV]Q.{3}[RK]G.{3}[RK]..[RILVWY]'
    matching protein names (as a set): {['SPO_RS04045', 'SPO_RS15950']}
3. C2H2 zinc-finger domain
    regular expression: r'C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H'
    matching protein names (as a set): {['SPO_RS21685']}
"""