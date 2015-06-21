#!/usr/bin/python

import sys
#from nltk.tree import *

# Tree manipulation

# Extract phrases from a parsed (chunked) tree
# Phrase = tag for the string phrase (sub-tree) to extract
# Returns: List of deep copies;  Recursive
def ExtractPhrases( myTree, phrase):
    myPhrases = []
    if (myTree.node == phrase):
        myPhrases.append( myTree.copy(True) )
    for child in myTree:
        if (type(child) is Tree):
            list_of_phrases = ExtractPhrases(child, phrase)
            if (len(list_of_phrases) > 0):
                myPhrases.extend(list_of_phrases)
    return myPhrases

#######################################



def main():
    f = open(sys.argv[1], 'r+');
    lines = f.readlines()

    for line in lines:
        test = Tree.parse(line)
        print "Input tree: ", test

        print "\nNoun phrases:"
        list_of_noun_phrases = ExtractPhrases(test, 'NP')
        for phrase in list_of_noun_phrases:
            print " ", phrase

if __name__ == "__main__":
    main()
