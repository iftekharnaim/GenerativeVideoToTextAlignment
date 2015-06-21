#!/usr/bin/python
#
# Author: Iftekhar Naim, CS Dept., University of Rochester.
# Email: inaim@cs.rochester.edu
#

"""Train a HMM model for aligning parallel video and text data.

To run this code, please use this command:

> python hmm_multivideos.py [FLAG_HMM2] [FLAG_ANVIL]

The two flags are optional, and takes boolean values (True/False).
If not provided, then the default values (True) are used.
"""
import sys
from collections import defaultdict
import math
import operator
from log_rescaling import *
from time import *

# load custom class
from video_protocol_data import VideoProtocolData
from timing_info import TimeInterval
from input_file_processor import *
from global_variables import *
from hmm_library import *
from scipy.special import digamma, betaln, gammaln
############################
def ReadAndSetFlags():
  global FLAG_HMM2
  global FLAG_ANVIL
  # Process the input whether to use HMM1 or HMM2. Default value is True.
  if len(sys.argv) > 1:
    if sys.argv[1] == "True":
      FLAG_HMM2 = True
    elif sys.argv[1] == "False":
      FLAG_HMM2 = False
    else:
      print("Incorrect argument: the first argument must be True/False")
      exit(-1)
  # ANVIL or Vision data? Default is True, i.e., ANVIL.
  if len(sys.argv) > 2:
    if sys.argv[2] == "True":
      FLAG_ANVIL = True
    elif sys.argv[2] == "False":
      FLAG_ANVIL = False
    else:
      print("Incorrect argument: the second argument must be True/False")
      exit(-1)
############################
def main():
  global next_states
  global prev_states
  global FLAG_HMM2
  global FLAG_ANVIL
  global FLAG_VB
  
  global protocols
  global anvil_tags
  global vision_tags

  ReadAndSetFlags()
  UpdateFlagsHMM(prev_states, next_states, FLAG_HMM2, FLAG_ANVIL, FLAG_VB)
  print FLAG_HMM2, FLAG_ANVIL
  # Read the input sequences in the data items
  protocols = protocols_anvil if FLAG_ANVIL else protocols_vision
  data_items = ReadDataItems(protocols, anvil_tags, vision_tags)
  Xg, Yg, Zg = GetGlobalCounts(data_items)
  
  print len(Xg), len(Yg), len(Zg)
  PrintGlobalDictionary(Xg, Yg, Zg, data_items)
  EM_START_TIME = time()
  # First Fit HMM1 model
  UpdateFlagsHMM(prev_states, next_states, FLAG_HMM2, FLAG_ANVIL, FLAG_VB)
  M, C = FitHMMAlignmentModel(data_items, Xg, Yg, nIterationsEM)
  ########## EM Iterations ended #############
  EM_ELAPSED_TIME = (time() - EM_START_TIME)
  PrintModelParametersHMM(M, C, Xg, Yg, 3)
  # print alignment accuracies
  accuracies, avg_acc = ViterbiAlignmentAccuraciesHMM(data_items, M, C)
  # print matching accuracies
  MatchingAccuraciesHMM(data_items, M, C)

  print("---------")
  print("FLAG HMM2: " + str(FLAG_HMM2))
  print("Elapsed time: " + str(EM_ELAPSED_TIME) + " secs")
  print "Average acc: %f" % avg_acc

if __name__ == "__main__":
  main()