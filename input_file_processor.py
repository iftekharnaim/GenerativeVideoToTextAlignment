#!/usr/bin/python
#
# Author: Iftekhar Naim, CS Dept., University of Rochester.
# Email: inaim@cs.rochester.edu
#

import sys
from nltk.tree import *
from collections import defaultdict
import math
import operator
from log_rescaling import *
# load custom class
from timing_info import TimeInterval

#######################################
# Extract phrases from a parsed (chunked) tree
# Phrase = tag for the string phrase (sub-tree) to extract
# Returns: List of deep copies;  Recursive
def ExtractPhrases(myTree, phrase):
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
# Extract phrases from a parsed (chunked) tree
# Phrase = tag for the string phrase (sub-tree) to extract
# Returns: List of deep copies;  Recursive
def ExtractNounPhrases(myTree):
  myPhrases = []
  if (myTree.node == "NN") | (myTree.node == "NNP"):
    myPhrases.append( myTree.copy(True) )
  
  prev_child = None
  next_child = None
  
  child_index = -1
  
  for child in myTree:
    child_index = child_index + 1
    
    if (type(child) is Tree):
      if (child_index < len(myTree) - 1):
        next_child = myTree[child_index + 1]
      
      if prev_child == None:
        if child_index >= (len(myTree)-1):
          list_of_phrases = ExtractNounPhrases(child)
        elif (next_child.node  == "NN") | (next_child.node == 'NNP') | ((next_child.node == "PP") & (next_child.leaves()[0].lower() == "of") ):
          list_of_phrases = []
        else:
          list_of_phrases = ExtractNounPhrases(child)
      elif (prev_child.node == "VB") & (len(prev_child.leaves()) == 1) & (prev_child.leaves()[0].lower() == 'write'):
        list_of_phrases = []
      else:
        if child_index >= (len(myTree)-1):
          list_of_phrases = ExtractNounPhrases(child)
        elif (next_child.node  == "NN") | (next_child.node == 'NNP') | ((next_child.node == "PP") & (next_child.leaves()[0].lower() == "of") ):
          list_of_phrases = []
        else:
          list_of_phrases = ExtractNounPhrases(child)
      
      prev_child = child
      
      if (len(list_of_phrases) > 0):
        myPhrases.extend(list_of_phrases)
  return myPhrases

#######################################
# Read the noun phrases from the parse tree
# Input: the input file with parse trees

def GetPrepositionPhrasesFromNotes(parse_file):
  # open the file containing all the parses
  f = open(parse_file, 'r+');
  lines = f.readlines()
  f_notes = open(parse_file.replace('parses','protocol'), 'r+')
  lines_notes = f_notes.readlines()
  
  preps = [];
  step_index = [];
  
  f.close()
  f_notes.close()
  
  count = -1;
  
  for line in lines:
    count = count + 1;
    # read the bracketed parse tree, and create a tree object
    tree = Tree.parse(line)
    
    preps_line = ExtractPhrases (tree, "PP")
    #nouns_in_this_line = ["null" + str(count/2)]
    preps_in_this_line = []
    for phrase in preps_line:
      phrase = str(phrase)
      phrase = phrase.replace(')', '')
      phrase_parts = phrase.split(" ")
      #            phrase = phrase_parts[1].lstrip().lower();
      
      # ignore when "step" is detected as noun
      if phrase.lower() == "step":
        continue;
      
      #            print phrase
      if phrase in preps_in_this_line:
        continue
      
      preps_in_this_line.append(phrase)
    #if(len(preps_in_this_line)>0):
    preps.append(preps_in_this_line);
  
  step_count = 0;
  
  for line in lines_notes:
    line = line.lower()
    line = line.lstrip().rstrip()
    parts = line.split(' ');
    if len(parts) <= 1:
      continue
    
    if parts[0].lower() == 'step':
      step_count = step_count + 1
    step_index.append(step_count)
  
  return (preps, lines,step_index);

#######################################
# Read the verbs from the parse tree
# Input: the input file with parse trees
def GetVerbsFromNotes(parse_file):
  # open the file containing all the parses
  f = open(parse_file, 'r+');
  lines = f.readlines()
  f_notes = open(parse_file.replace('parses','protocol'), 'r+')
  lines_notes = f_notes.readlines()
  f.close()
  f_notes.close()
  
  verbs = [];
  step_index = [];
  count = -1;
  
  for line in lines:
    count = count + 1;
    # read the bracketed parse tree, and create a tree object
    tree = Tree.parse(line)
    
    verbs_line = ExtractPhrases (tree, "VB")
    verbs_in_this_line = []
    for phrase in verbs_line:
      phrase = str(phrase)
      phrase = phrase.replace(')', '')
      phrase_parts = phrase.split(" ")
      phrase = phrase_parts[1].lstrip().lower();
      if phrase in verbs_in_this_line:
        continue
      verbs_in_this_line.append(phrase)
    verbs.append(verbs_in_this_line);
  
  step_count = 0;
  for line in lines_notes:
    line = line.lower()
    line = line.lstrip().rstrip()
    parts = line.split(' ');
    if len(parts) <= 1:
      continue
    
    if parts[0].lower() == 'step':
      step_count = step_count + 1
    step_index.append(step_count)
  
  return (verbs, lines, step_index);
#######################################
# Read the noun phrases from the parse tree
# Input: the input file with parse trees
def GetNounPhrasesFromNotes(parse_file):
  # open the file containing all the parses
  f = open(parse_file, 'r+');
  lines = f.readlines()
  f_notes = open(parse_file.replace('parses','protocol'), 'r+')
  lines_notes = f_notes.readlines()
  f.close()
  f_notes.close()
  
  nouns = []
  step_index = []
  count = -1
  
  for line in lines:
    count = count + 1;
    # read the bracketed parse tree, and create a tree object
    tree = Tree.parse(line)
    nouns_line = ExtractNounPhrases(tree)
    nouns_in_this_line = []
    for phrase in nouns_line:
      phrase = str(phrase)
      phrase = phrase.replace(')', '')
      phrase_parts = phrase.split(" ")
      phrase = phrase_parts[1].lstrip().lower();
      
      # ignore when "step" is detected as noun
      if phrase == "step":
        continue;
      
      if phrase in nouns_in_this_line:
        continue
      
      nouns_in_this_line.append(phrase)
    if(len(nouns_in_this_line)>0):
      nouns.append(nouns_in_this_line);
  
  step_count = 0;
  
  for line in lines_notes:
    line = line.lower()
    line = line.lstrip().rstrip()
    parts = line.split(' ');
    if len(parts) <= 1:
      continue
    if parts[0].lower() == 'step':
      step_count = step_count + 1
    step_index.append(step_count)
  
  return (nouns, lines,step_index);


###################################
# Read the Anvil tags from the file
def GetAnvilTags(anvil_file):
  """Reads the blob tracks from the input file."""
  f = open(anvil_file, 'r')
  lines = f.readlines();
  f.close();
  
  anvilTags = [];
  lines_avl = [];
  time_intervals = [];
  count = -1;
  step_index_avl = [];
  
  for line in lines:
    count = count + 1;
    if ( count == 0):
      continue;
    
    # skip the intervals without any anvil tags
    parts = line.split(" ");
    if len(parts) < 4:
      continue;
    tags = parts[3].lstrip().rstrip();
    if tags == "":
      continue;
    
    # get the starting time of the interval
    start_time = float(parts[1])
    end_time = float(parts[2])
    mean_time = start_time + (end_time - start_time)/2.0;
    duration = end_time - start_time
    
    interval = 1.0;
    if(duration > interval) :
      val = duration/interval + 1.0;
      count = int(val)
      for window in range(0,count):
        tags_list = tags.split(",");
        anvilTags.append(tags_list[:-1])
        lines_avl.append(line);
        step_index_avl.append(int(tags_list[-1]))
        time_intervals.append( TimeInterval(start_time + window * interval, start_time + (window+1) * interval, start_time + (window+0.5) *interval));
    else:
      tags_list = tags.split(",");
      step_index_avl.append(int(tags_list[-1]))
      anvilTags.append(tags_list[:-1])
      lines_avl.append(line);
      time_intervals.append( TimeInterval(start_time, end_time, mean_time));

  return (anvilTags, lines_avl, time_intervals, step_index_avl);

################################
def GetCountAnvilTags(tags, time_info_anvil):
  """Returns a dictionary of the counts of the video blobs/tags."""
  Y = {}
  timeY = {}
  count = 0
  for taglist in tags:
    for tag in taglist:
      Y[tag] = Y.setdefault(tag, 0) + 1
      timeY.setdefault(tag, []).append(time_info_anvil[count])
    
    count = count + 1
  return (Y, timeY)

################################
def GetCountNounTags(nouns, time_info_nouns):
  """Returns a dictionary of the counts of the nouns."""
  # set of distinct objects
  X = {}
  timeX = {}
  count = 0
  for nounlist in nouns:
    for noun in nounlist:
      X[noun] = X.setdefault(noun, 0) + 1
      timeX.setdefault(noun, []).append(time_info_nouns[count])
    count = count + 1
  return X, timeX

################################
def GetCountVerbTags(verbs):
  """Returns a dictionary of counts of the verbs."""
  # set of distinct objects
  Z = {}
  for verblist in verbs:
    for verb in verblist:
      Z[verb] = Z.setdefault(verb, 0) + 1

  return Z

############################
def GetNounPhrasesAndVerbsFromNotes(parse_file):
  """Extracts verbs from the protocol parse."""
  (verbs, lines_text,step_index_nouns) = GetVerbsFromNotes(parse_file);
  (nouns, lines_text,step_index_nouns) = GetNounPhrasesFromNotes(parse_file);
  for i in range(len(verbs)):
    for j in range(len(verbs[i])):
      if verbs[i][j] not in nouns[i]:
        nouns[i].append(verbs[i][j])
  return nouns, lines_text, step_index_nouns

########### helper functions ##########
def PrintParseTrees(file):
  """Print all the parse trees in the input parse file."""
  f = open(file, 'r')
  lines = f.readlines()
  f.close()
  for line in lines:
    test = Tree.parse(line)
    print "Input tree: ", test

############################
# Acknowledgement: The code for powerset generation is copied from:
# http://www.technomancy.org/python/powerset-generator-python/
def powerset(seq):
  """Returns all the subsets of a set."""
  if len(seq) <= 1:
    yield seq
    yield []
  else:
    for item in powerset(seq[1:]):
      yield [seq[0]]+item
      yield item
