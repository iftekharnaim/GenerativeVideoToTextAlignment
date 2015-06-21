#!/usr/bin/python
#
# Author: Iftekhar Naim, CS Dept., University of Rochester.
# Email: inaim@cs.rochester.edu
#

"""Train an HMM model with latent observation parameters, for aligning parallel video and text data.
  
  To run this code, please use this command:
  
  > python latent_hmm_multivideos.py [FLAG_HMM2] [FLAG_ANVIL]
  
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

############################################
#def DecodingGivenMatching(nouns, tags, X, Y, protocol):
#  M = defaultdict(dict)
#  
#  for x in X:
#    for y in Y:
#      M[x][y] = 0.0
#  
#  if (protocol == 'YPAD'):
#    # load the YPAD ground truth assignments
#    print('YPAD')
#  elif (protocol == 'LLGM'):
#    # load the LLGM assignemnt
#    print('LLGM')
#  else:
#    # load the CELLOBIOSE assignment
#    print('')
#

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
  PrintGlobalDictionary(Xg, Yg, Zg, data_items)
  # First initialize using HMM1 model
  nIterations = 15 if FLAG_HMM2 else 50
  next_states = range(0,2,1)
  prev_states = range(-1,1,1)
  UpdateFlagsHMM(prev_states, next_states, FLAG_HMM2, FLAG_ANVIL, FLAG_VB)
  M, C = InitializeUniform(Xg, Yg)
  M, C = ApplyEM(data_items, M, C, Xg, Yg, nIterations)
  # Fit HMM2, if FLAG_HMM2 is true
  if FLAG_HMM2:
    print "Start HMM2 iterations"
    nIterations = 35
    next_states = range(-2,3,1)
    prev_states = range(-2,3,1)
    UpdateFlagsHMM(prev_states, next_states, FLAG_HMM2, FLAG_ANVIL, FLAG_VB)
    # reinitialize the transition probabilities
    C = InitializeUniformTransitionProbs()
    M, C = ApplyEM(data_items, M, C, Xg, Yg, nIterations)

  # initialize observation parameters
  O = UniformObservationParameters(Xg, 0.9)
  EM_START_TIME = time()
  M, C, O = ApplyLatentObservationEM(data_items, M, C, O, Xg, Yg, nIterationsEM)
  EM_ELAPSED_TIME = (time() - EM_START_TIME)
  ########## EM Iterations ended #############
  PrintModelParametersHMM(M, C, Xg, Yg, 3)
  # needs to be fixed
  accuracies, avg_acc = ViterbiAlignmentAccuraciesLatentHMM(data_items, M, C, O)
  print "-----------"
  for x in Xg:
    print(x,O[x])

  print "---------"
  print "FLAG HMM2: " + str(FLAG_HMM2)
  print "Elapsed time: " + str(EM_ELAPSED_TIME) + " secs"
  print "Time per iteration per video:" + str(float(EM_ELAPSED_TIME)/(float(nIterationsEM)*float(len(data_items))))
  print "Average acc: %f" % avg_acc

if __name__ == "__main__":
  main()
