#!/usr/bin/python
#
# Author: Iftekhar Naim, CS Dept., University of Rochester.
# Email: inaim@cs.rochester.edu
#

"""Library functions for HMM/IBM1 model."""
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
from scipy.special import digamma, betaln, gammaln

################################
def IsValidState(i,j, m, n):
  """Returns true if the state (i,j) is a valid state.
    
  Inputs: 
    A state (i,j) and the lattice dimension (m,n).
  Output:
    Returns True if (i,j) is a valid state, given the dimension.
  """
  if (i > m) or (j > n) or (i <= 0) or (j <= 0):
    return False
  return True
################################
def GetBackPointers(i, j, m, n):
  """Returns the list of backpointers for a given state (i,j)."""
  global prev_states
  BP = list()
  if j == 1:
    b = (i-1,0)
    BP.append(b)
    return BP
  
  for di in prev_states:
    if IsValidState(i+di,j-1, m, n):
      BP.append((i+di,j-1))
  return BP
################################
def GetSuccessors(i, j, m, n):
  """Returns the successors states of (i,j)."""
  global next_states
  Sc = list()
  for di in next_states:
    if IsValidState(i+di, j+1, m, n):
      Sc.append((i+di,j+1))
  return Sc
################################
def PrintBestAlignment(nouns, tags, B,D):
  """Prints the best alignment, given the viterbi matrices."""
  m = len(nouns)
  n = len(tags)
  print("Print the best alignment")
  
  end_point = -1
  max_end = -1e100
  
  for k in range(1,m+1):
    if(D[k][n] > max_end):
      max_end = D[k][n]
      end_point = k
  
  (i,j) = (k,n)
  
  BestAlignment = []
  BestAlignment.append((i,j))
  
  output_text = "";
  str1 = "\n==========\n";
  str1 = str1 + str(i-1) + ": Nouns in sentence: " + str(nouns[i-1]) + "\n";
  str1 = str1 + str(j-1) + ": Anvil tags: " +  str(tags[j-1]) + "\n";
  output_text = str1;
  
  while(True):
    (i,j) = B[i][j]
    if (j<=0):
      break;
    BestAlignment.append((i,j))
    str1 = ""
    str1 = str1 + str(i-1) + ": Nouns in sentence: " + str(nouns[i-1]) + "\n";
    str1 = str1 + str(j-1) + ": Anvil tags: " +  str(tags[j-1]) + "\n";
    str1 = str1 + "==========\n";
    output_text = str1 + output_text;
  
  output_text = "\n\n\n\n" + output_text
  print(output_text)
  return (BestAlignment,output_text)
################################
def GetBestAlignment(nouns, tags, B,D):
  """Returns the best alignment, given the viterbi matrices.
     It is similar to the method PrintBestAlignment, except it only returns
     the best alignment, and does not print it to stdout.
  """
  m = len(nouns)
  n = len(tags)
  
  end_point = -1
  max_end = -1e200
  
  for k in range(1,m+1):
    if(D[k][n] > max_end):
      max_end = D[k][n]
      end_point = k
  
  (i,j) = (k,n)
  
  BestAlignment = []
  BestAlignment.append(i-1)
  
  while(True):
    (i,j) = B[i][j]
    if (j<=0):
      break;
    BestAlignment.append(i-1)
  BestAlignment.reverse()
  #print("Return from the best alignment get method")
  #print(BestAlignment)
  return BestAlignment

#######################################
def IBM1Prob(nouns, tags, M):
  """Returns the IBM Model 1 probability for a given set of nouns and blobs."""
  p = 1.0;
  for tag in tags:
    sum_tag = 0.0
    for noun in nouns:
      sum_tag = sum_tag + M[noun][tag]
    p = p * sum_tag;
  mm = len(tags)
  ll = len(nouns)
  scale = ll ** mm;
  return p/float(scale);
#######################################
def GetIBM1Probabilities(nouns, tags, M):
  """Get the IBM1 probabilities for all the pairwise alignments."""
  # f = anvil tag (i.e. blob)
  # e = noun
  m = len(nouns)
  n = len(tags)
  S = [[0 for j in range(n+1)] for i in range(m+1)]
  # Initialize the viterbi scores
  for i in range(1,m+1):
    S[i][0] = 0.0
  for j in range(1,n+1):
    S[0][j] = 0.0
  
  # Viterbi dynamic programming
  for i in range(1,m+1):
    for j in range(1,n+1):
      S[i][j] = IBM1Prob(nouns[i-1],tags[j-1],M);
  return S
#######################################
def NormalizedExpectedCounts(nouns, tags, M):
  """Normalize the expected counts for the E-step.
    This is necessary because each tag/blob came from one of the nouns. So the
    probability should add up to 1 when summed over the nouns.
    The function is called from the method GetExpectedCounts().
  """
  EC_ij = defaultdict(dict)
  
  for nn in nouns:
    for tt in tags:
      EC_ij[nn][tt] = 0.0
  
  p = 1.0;
  for tag in tags:
    sum_tag = 0.0
    for noun in nouns:
      sum_tag = sum_tag + M[noun][tag]
    for noun in nouns:
      EC_ij[noun][tag] = M[noun][tag]/sum_tag
  return EC_ij;

#######################################
def GetExpectedCounts(nouns, tags, M):
  """Returns the expected counts for all pairwise alignments."""
  # f = tag
  # e = noun
  m = len(nouns)
  n = len(tags)

  EC = [[defaultdict(dict) for j in range(n+1)] for i in range(m+1)]
  # Initialize the expected count dict in each position
  for i in range(1,m+1):
    EC[i][0] = defaultdict(dict)
  for j in range(1,n+1):
    EC[0][j] = defaultdict(dict)
  # estimate the expected counts
  for i in range(1,m+1):
    for j in range(1,n+1):
      EC[i][j] = NormalizedExpectedCounts(nouns[i-1],tags[j-1],M)
  return EC
################################
def forward_backward(nouns, tags, M, C):
  """The method for estimating forward and backward probabilities and
    the expected counts for the E-step.
    
    This is one of the most important methods in this project.
  """
  m = len(nouns)
  n = len(tags)
  global next_states
  global prev_states
  global FLAG_HMM2
  # Estimate IBM1 probabilities for each pairwise alignment
  S = GetIBM1Probabilities(nouns, tags, M)
  # Estimate the expected counts
  EC = GetExpectedCounts(nouns, tags, M)
  #######################################
  # Perform the forward loop
  FW = [[eln(0.0) for j in range(n+1)] for i in range(m+1)]
  # Initialize the forward scores
  FW[0][0] = eln(1.0);

  # Forward recursion
  for j in range(1,n+1):
    for i in range(1,m+1):
      fw_sum = eln(0.0)
      BP = GetBackPointers(i, j, m, n)
      for k in range(0, len(BP)):
        i1 = BP[k][0]
        j1 = BP[k][1]
        jump = i - i1
        fw_sum = elnsum(fw_sum, elnproduct(FW[i1][j1], eln(C[str(jump)])))
      fw_sum = elnproduct(fw_sum, eln(S[i][j]))
      FW[i][j] = fw_sum
  #######################################
  # Perform backward loop
  BW = [[eln(0.0) for j in range(n+1)] for i in range(m+1)]
  for i in range(m,-1,-1):
#    BW[i][n] = eln(1.0)
    jump_key = str(m-i)
    if jump_key in C:
      BW[i][n] = eln(1.0)
    else:
      BW[i][n] = eln(0.0)
  for j in range(n-1,-1,-1):
    for i in range(m,-1,-1):
      bw_sum = eln(0.0);
      Sc = GetSuccessors(i, j, m, n);
      for k in range(0, len(Sc)):
        i1 = Sc[k][0]
        j1 = Sc[k][1]
        jump = i1 - i
        bw_sum = elnsum(bw_sum, elnproduct(eln(C[str(jump)]),elnproduct(BW[i1][j1], eln(S[i1][j1]))) );
      BW[i][j] = bw_sum
  #######################################
  # Estimate the gamma values
  Gamma = [[eln(0.0) for j in range(n+1)] for i in range(m+1)]
  for j in range(0,n+1):
    Z = eln(0.0)
    for i in range(0,m+1):
      Gamma[i][j] = elnproduct(FW[i][j], BW[i][j])
      Z = elnsum(Z, Gamma[i][j])
    Zinv = Z * -1.0
    for i in range(0,m+1):
      Gamma[i][j] = eexp(elnproduct(Gamma[i][j], Zinv));
  ##########################################################
  # Estimate the expected counts for jump probabilities
  C_item = dict(C)
  # initialize expected counts for jumps to zero
  for i in next_states:
    C_item[str(i)] = eln(0.0)
  for j in range(0,n+1):
    tempZ = eln(0.0)
    temp_C = {}
    for i in range(0,m+1):
      Sc = GetSuccessors(i, j, m, n)
      for k in range(len(Sc)):
        i1 = Sc[k][0]
        j1 = Sc[k][1]
        jump = i1 - i
        key = str(jump)
        val = elnproduct(elnproduct(eln(S[i1][j1]), eln(C[key])), elnproduct(FW[i][j], BW[i1][j1]));
        if key in temp_C:
          temp_C[key] = elnsum(temp_C[key],val)
        else:
          temp_C[key] = val
        tempZ = elnsum(tempZ , val)
    # after the loop over i.
    # Normalize the expected counts for the current j
    tempZ = tempZ * -1.0
    for key in temp_C:
      temp_C[key] = elnproduct(temp_C[key] , tempZ)
      # add up the expected counts for the jump size key
      C_item[key] = elnsum(C_item[key], temp_C[key])

  return (FW,BW,Gamma, C_item, EC)
################################
def FitHMMAlignmentModel(data_items, Xg, Yg, nIterations):
  global next_states
  global prev_states
  # start with HMM1
  next_states = range(0,2,1)
  prev_states = range(-1,1,1)
  M, C = InitializeUniform(Xg, Yg)
  if not FLAG_HMM2:
    M, C = ApplyEM(data_items, M, C, Xg, Yg, nIterations)
  else:
    M, C = ApplyEM(data_items, M, C, Xg, Yg, 10)
    # now run HMM2 model fitting
    next_states = range(-2,3,1)
    prev_states = range(-2,3,1)
    # reinitialize the transition probabilities
    C = InitializeUniformTransitionProbs()
    M, C = ApplyEM(data_items, M, C, Xg, Yg, nIterations-10)
  return M, C

################################
def ApplyEM(data_items, M, C, Xg, Yg, nIterations):
  for iter in range(nIterations):
    ECList = []
    GammaList = []
    CList = []
    print("Iteration: " + str(iter))
    # E-step, by iterating over all the data items
    for item in data_items:
      FW, BW, Gamma, C_item, EC = forward_backward(item.nouns, item.tags, M, C)
      GammaList.append(Gamma)
      ECList.append(EC)
      CList.append(C_item)
    # M-step
    M, C = Mstep(data_items, ECList, CList, GammaList, Xg, Yg)
  return M, C
################################
def InitializeUniform(X, Y):
  """Initializes the parameters M and C uniformly."""
  global next_states
  # Initialize cooccurrence probabilities
  M = defaultdict(dict)
  for x in X:
    for y in Y:
      M[x][y] = 1.0 / float(len(Y))
  # Initialize the transition parameters
  C = {}
  for i in next_states:
    C[str(i)] = (1 / float(len(next_states)))
  return M, C

################################
def InitializeUniformTransitionProbs():
  """Uniformly initializes the transition probability."""
  global next_states
  # The transition parameters
  C = {}
  for i in next_states:
    C[str(i)] = 1.0 / float(len(next_states))
  return C
################################
def Mstep(data_items, ECList, CList, GammaList, Xg, Yg):
  """Normalize the expected counts for the M-step."""
  global FLAG_VB
  global next_states
  global prev_states
  
  M = defaultdict(dict);
  # initialize to zero
  for x in Xg:
    for y in Yg:
      M[x][y] = 0.0;

  for index in range(len(data_items)):
    nouns = data_items[index].nouns
    tags = data_items[index].tags
    Gamma = GammaList[index]
    EC = ECList[index]
    m = len(nouns)
    n = len(tags)
    for i in range(1,m+1):
      for j in range(1,n+1):
        nn = nouns[i-1];
        tt = tags[j-1];
        for noun in nn:
          for tag in tt:
            M[noun][tag] = M[noun][tag] + (Gamma[i][j] * EC[i][j][noun][tag])

  dir_alpha = 1.0
  if FLAG_VB == False:
    dir_alpha = 0.0
  K = float(len(Yg))
  K_alpha = K * dir_alpha
  
  for x in Xg:
    Z = 0.0
    for y in Yg:
      Z = Z + M[x][y];
    
    for y in Yg:
      if FLAG_VB == True:
        if(math.exp(digamma(Z + K_alpha)) < 1e-40):
          M[x][y] = 1e-40
        else:
          M[x][y] = math.exp(digamma( M[x][y]+dir_alpha) )  / math.exp(digamma(Z + K_alpha));
      else:
        M[x][y] = M[x][y]/ Z;
      if M[x][y] < 1e-40:
        M[x][y] = 1e-40

  # udpate the C probabilities
  # Normalize the jump probabilities, this should be moved to the M-step part
  Z_c = eln(0.0)
  # initialize to zero
  C = {}
  for i in next_states:
    C[str(i)] = eln(0.0)
  for index in range(len(data_items)):
    C_item = CList[index]
    for i in next_states:
      key = str(i)
      C[key] = elnsum(C[key], C_item[key])
      Z_c = elnsum(Z_c, C_item[key])
  Z_c_inv = Z_c * -1.0

  for i in next_states:
    key = str(i)
    C[key] = eexp(elnproduct(C[key], Z_c_inv))
    if (C[key] < 1e-40 ):
      C[key] = 1e-40

  print("Current C value:")
  print(C)
  return M, C

###########################################
def JointViterbi(nouns, tags, X, Y, M, C):
  # perform a joint decoding over all the states and their associated IBM-1 probs
  m = len(nouns)
  n = len(tags)
  # array for tracking the best path in dynamic programming
  D = [[1e-200 for j in range(n+1)] for i in range(m+1)]
  # the array of best backpointers for each state
  B= [[(-1,-1) for j in range(n+1)] for i in range(m+1)]

  S = GetIBM1Probabilities(nouns, tags, M);
  Match = defaultdict(dict)
  
  for x in X:
    for y in Y:
      Match[x][y] = 0.0
  
  # Initialize the viterbi scores
  for i in range(1,m+1):
    D[i][0] = -1e200;
  
  for j in range(1,n+1):
    D[0][j] = -1e200;
  
  # Viterbi dynamic programming
  for j in range(1,n+1):
    for i in range(1,m+1):
      nn = nouns[i-1];
      tt = tags[j-1];
      
      BP = GetBackPointers(i,j, m, n)
      max_val = -1e2000
      jump_size = -1000
      for k in range(len(BP)):
        i1 = BP[k][0]
        j1 = BP[k][1]
        
        jump_size = i - i1
        temp_val = D[i1][j1] + math.log(S[i][j]+1e-200) + math.log(C[str(jump_size)]+1e-200)
        if  temp_val >= max_val:
          max_val = temp_val
          B[i][j] = (i1,j1)
          jump_size = i-i1
      
      D[i][j]  = max_val
  
  return (D,B, Match)

###########################################
def UniformBaseline(timing_info_anvil, nouns, tags, step_index_nouns, step_index_avl):
  i = 0
  j = 0
  
  total_interval = timing_info_anvil[-1].end - timing_info_anvil[0].start;
  offset = timing_info_anvil[i].start;
  
  time_per_sentence = total_interval / float(len(nouns));
  tags_per_sentence = int(float(len(tags)) / float(len(nouns)))+1;
  
  
  alignment = []
  for i in range(len(tags)):
    tag_assignment = int(float(i) / float(tags_per_sentence)) + 1
    alignment.append((i+1,tag_assignment))
  
  # first read the ANVIL file
  count = 0
  error_count = 0
  
  for (i,j) in alignment:
    count = count + 1
    if step_index_avl[i-1] == step_index_nouns[j-1]:
      continue
    error_count = error_count + 1
  
  print("Baseline----");
  print(1.0-float(error_count)/float(count))

####################################
def PrintGlobalDictionary(Xg, Yg, Zg, data_items):
  print "\nNoun dictionary ---"
  for x in sorted(Xg):
    xstring = [x, ":", str(Xg[x])]
    for data in data_items:
      if x in sorted(data.X):
        value = data.X[x]
      else:
        value = 0
      xstring.append("|%d" % value)
    print(" ".join(xstring))

  print "\nBlob dictionary ---"
  for y in sorted(Yg):
    ystring = [y, ":", str(Yg[y])]
    for data in data_items:
      if y in data.Y:
        value = data.Y[y]
      else:
        value = 0
      ystring.append("|%d" % value)
    print(" ".join(ystring))

  print "\nVerb dictionary ---"
  for z in sorted(Zg):
    zstring = [z, ":", str(Zg[z])]
    for data in data_items:
      if z in data.Z:
        value = data.Z[z]
      else:
        value = 0
      zstring.append("|%d" % value)
    print(" ".join(zstring))

##################################
def GetTimingInfoNounsUniform(timing_info_anvil, num_sentences):
  # get the timing info for the text nouns
  timing_info_nouns = []
  count = 0.0;
  video_duration = timing_info_anvil[-1].end
  for i in range(num_sentences):
    count = count + 1.0;
    tval = video_duration * count / float(num_sentences)
    timing_info_nouns.append(TimeInterval(tval,tval,tval))
  return timing_info_nouns
###############################
def GetGlobalCounts(data_items):
  """Combines the per-video dictionaries to get global dictionaries and counts."""
  Xg = {}  # noun dictionary
  Yg = {}  # blob dictionary
  Zg = {}  # verb dictionary
  for data in data_items:
    for x in data.X:
      Xg[x] = Xg.setdefault(x, 0) + data.X[x]
    for y in data.Y:
      Yg[y] = Yg.setdefault(y, 0) + data.Y[y]
    for z in data.Z:
      Zg[z] = Zg.setdefault(z, 0) + data.Z[z]

  return Xg, Yg, Zg
#############################
def StepwiseAlignmentAccuracy(alignment, step_index_avl, step_index_nouns):
  count = 0
  error_count = 0
  for i in range(len(alignment)):
    count = count + 1
    j = alignment[i]
    if step_index_avl[i] == step_index_nouns[j]:
      continue
    error_count = error_count + 1
  return 1.0-float(error_count)/float(count)

############################
def PrintBestMatching(nouns, tags, bestAlignment, M, X, Y):
  for y in Y:
    best_p = 0
    best_j = 0
    for x in X:
      prob = 0.0
      # iterate over all the alignments
      for i in range(len(bestAlignment)):
        j = bestAlignment[i]
        if (y not in tags[i]) | (x not in nouns[j]):
          continue
        
        prob = prob + M[x][y]
        if prob >= best_p:
          best_p = prob;
          best_j = x
    print(y,best_j)

###############################
def ReadDataItems(protocols, anvil_tags, vision_tags):
  """Read protocol and video data from the input files."""
  global FLAG_ANVIL
  data_items = []
  
  # do some sanity checks
  if FLAG_ANVIL:
    assert len(protocols) == len(anvil_tags)
  else:
    assert len(protocols) == len(vision_tags)

  # iterate over all the protocol entries and read data
  for data_index in range(len(protocols)):
    parse_file = "./protocols/parses_%s.txt" % protocols[data_index]
    if FLAG_ANVIL == True:
      anvil_file = "./video_blobs/%s" % anvil_tags[data_index]
    else:
      anvil_file = "./video_blobs/%s" % vision_tags[data_index]
    
    nouns, lines_text, step_index_nouns = GetNounPhrasesFromNotes(parse_file)
    verbs, lines_text, step_index_dummy = GetVerbsFromNotes(parse_file)
    tags, lines_avl, timing_info_anvil,step_index_avl = GetAnvilTags(anvil_file)
    
    timing_info_nouns = GetTimingInfoNounsUniform(timing_info_anvil, len(nouns))
    X, timeX = GetCountNounTags(nouns, timing_info_nouns)
    Y, timeY = GetCountAnvilTags(tags, timing_info_anvil)
    Z = GetCountVerbTags(verbs)
    
    dataitem = VideoProtocolData(nouns, tags, verbs, step_index_nouns, step_index_avl,lines_text, lines_avl, X, Y, Z, timeX, timeY, timing_info_anvil)
    data_items.append(dataitem)
    print(len(dataitem.nouns), len(dataitem.tags))
  return data_items
############################
def ViterbiAlignmentAccuraciesHMM(data_items, M, C):
  accuracies = []
  avg_acc = 0.0
  nItems = len(data_items)
  for item in data_items:
    D, B, Match = JointViterbi(item.nouns, item.tags, item.X, item.Y, M, C)
    bestAlignment = GetBestAlignment(item.nouns, item.tags, B, D)
    acc = StepwiseAlignmentAccuracy(bestAlignment, item.step_index_avl, item.step_index_nouns)
    print("Accuracy: %f" % acc)
    accuracies.append(acc)
    avg_acc = avg_acc + acc/float(nItems)
  print "Avg accuracy HMM: %f" % avg_acc
  return accuracies, avg_acc

############################
def MatchingAccuraciesHMM(data_items, M, C):
  """Prints the matching accuracies using IBM1/Viterbi Decoding"""
  for item in data_items:
    D, B, Match = JointViterbi(item.nouns, item.tags, item.X, item.Y, M, C)
    bestAlignment = GetBestAlignment(item.nouns, item.tags, B, D)
    PrintBestMatching(item.nouns, item.tags, bestAlignment, M, item.X, item.Y)
    print("==========")

############################
def PrintModelParametersHMM(M, C, Xg, Yg, K):
  """Prints the learned param values for HMM/IBM1."""
  print "Learned matching probs:"
  for x in sorted(Xg):
    print("========" + x + "==========")
    sortedX = sorted(M[x].iteritems(), key=operator.itemgetter(1), reverse=True)
    for i in range(K):
      print sortedX[i]

  print "Learned transition probs:"
  print(C)

############################
def UpdateFlagsHMM(pstates, nstates, fg_HMM2, fg_ANVIL, fg_VB):
  global prev_states
  global next_states
  global FLAG_VB
  global FLAG_HMM2
  global FLAG_ANVIL
  prev_states = pstates
  next_states = nstates
  FLAG_VB = fg_VB
  FLAG_HMM2 = fg_HMM2
  FLAG_ANVIL = fg_ANVIL
#################################################################
############################ Methods for Latent HMM #############
#################################################################
def KeyForNounSet(noun_set):
  """Creates a hashkey to identify the set of nouns used.
    
    The input is a set of nouns.
    Returns a string containing all the nouns, comma separated.
    """
  if len(noun_set) == 0:
    return ""
  key = noun_set[0]
  for i in range(1,len(noun_set)):
    key = key + "," + noun_set[i]
  return key

##############################
def ListOfNounsFromKey(key):
  """Returns a list of nouns from a given key (string of comma separated nouns)"""
  noun_set = []
  noun_set = key.split(",")
  return noun_set
#######################################
def ObservationProbabilityNounSubset(noun_list, noun_set, O):
  """Returns the observation probability.
    
  Input:
    A set of nouns, a subset of that set, and observation probabilities.
  Output:
    Returns the probability of observing that subset.
  """
  p = 1.0
  # multiply with the probability of the subset
  for nn in noun_list:
    if nn in noun_set:
      p = p * O[nn][1]
    else:
      p = p * O[nn][0]
  return p;
#######################################
def IBMProbabilitiesObservation(nouns, tags, M, O):
  m = len(nouns)
  n = len(tags)
  S = [[0 for j in range(n+1)] for i in range(m+1)]
  SP = [[{} for j in range(n+1)] for i in range(m+1)]
  # At each (i,j) we have a dictionary indexed by nounsets
  # At each of these dictionary entries, we'll have another 2D dictionary storing all the noun/tag prob pairs
  EC = [[{} for j in range(n+1)] for i in range(m+1)]
  
  # Initialize the viterbi scores
  for i in range(1,m+1):
    S[i][0] = 0.0
  for j in range(1,n+1):
    S[0][j] = 0.0

  # estimate the probabilities for each (i,j)
  for i in range(1,m+1):
    for j in range(1,n+1):
      noun_list = sorted(nouns[i-1])
      S[i][j] = 0.0
      # iterate over all the subsets in power set
      for noun_set in powerset(noun_list):
        # ignore the empty subset
        if len(noun_set) == 0:
          continue
        # get the IBM1 probablity for the given subset
        p = IBM1Prob(noun_set,tags[j-1],M) * ObservationProbabilityNounSubset(
            noun_list, noun_set, O)
        # get the noun key
        key = KeyForNounSet(noun_set)
        SP[i][j][key] = p
        S[i][j] = S[i][j] + p
        EC[i][j][key] = NormalizedExpectedCounts(noun_set,tags[j-1],M)

  return S, SP, EC

############################
def UniformObservationParameters(Xg, val):
  """Initialize noun observation parameters using val."""
  O = {}
  for noun in Xg:
    O[noun] = (1.0- val,val)
  return O
#########################################
def forward_backward_observation(nouns, tags, M, C, O):
  """Forward-backward E-step for latent HMM."""
  m = len(nouns)
  n = len(tags)
  global next_states
  global prev_states
  global FLAG_HMM2
  
  S, SP, EC = IBMProbabilitiesObservation(nouns, tags, M, O)
  ##### Perform the forward loop #####
  FW = [[eln(0.0) for j in range(n+1)] for i in range(m+1)]
  FW[0][0] = eln(1.0);
  # Forward recursion
  for j in range(1,n+1):
    for i in range(1,m+1):
      fw_sum = eln(0.0);
      BP = GetBackPointers(i, j, m, n)
      for k in range(0, len(BP)):
        i1 = BP[k][0]
        j1 = BP[k][1]
        jump = i - i1;
        fw_sum = elnsum(fw_sum,
                        elnproduct(FW[i1][j1], eln(C[str(jump)])))
      fw_sum = elnproduct(fw_sum, eln(S[i][j]))
      FW[i][j] = fw_sum
  #######################################
  # Perform backward loop
  BW = [[eln(0.0) for j in range(n+1)] for i in range(m+1)]
  for i in range(m,-1,-1):
#    BW[i][n] = eln(1.0)
    jump_key = str(m-i)
    if jump_key in C:
      BW[i][n] = eln(1.0)
    else:
      BW[i][n] = eln(0.0)

  for j in range(n-1,-1,-1):
    for i in range(m,-1,-1):
      bw_sum = eln(0.0);
      Sc = GetSuccessors(i, j, m, n)
      for k in range(0, len(Sc)):
        i1 = Sc[k][0]
        j1 = Sc[k][1]
        jump = i1 - i
        bw_sum = elnsum(bw_sum, elnproduct(eln(C[str(jump)]), elnproduct(BW[i1][j1], eln(S[i1][j1])) ) );
      BW[i][j] = bw_sum
  #######################################
  # Estimate the gamma values
  Gamma = [[eln(0.0) for j in range(n+1)] for i in range(m+1)]
  for j in range(0,n+1):
    Z = eln(0.0)
    for i in range(0,m+1):
      Gamma[i][j] = elnproduct(FW[i][j], BW[i][j])
      Z = elnsum(Z, Gamma[i][j])
    Zinv = Z * -1.0
    for i in range(0,m+1):
      Gamma[i][j] = eexp(elnproduct(Gamma[i][j], Zinv))
  ########################################
  # update the C values
  C_item = dict(C)
  # initialize expected counts for jumps to zero
  for i in next_states:
    C_item[str(i)] = eln(0.0)
  for j in range(0,n+1):
    tempZ = eln(0.0)
    temp_C = {}
    for i in range(0,m+1):
      Sc = GetSuccessors(i, j, m, n)
      for k in range(len(Sc)):
        i1 = Sc[k][0]
        j1 = Sc[k][1]
        jump = i1 - i
        key = str(jump)
        val = elnproduct(elnproduct(eln(S[i1][j1]), eln(C[key])), elnproduct(FW[i][j], BW[i1][j1]));
        if key in temp_C:
          temp_C[key] = elnsum(temp_C[key],val)
        else:
          temp_C[key] = val
        tempZ = elnsum(tempZ , val)
    # after the loop over i.
    # Normalize the expected counts for the current j
    tempZ = tempZ * -1.0
    for key in temp_C:
      temp_C[key] = elnproduct(temp_C[key] , tempZ)
      # add up the expected counts for the jump size key
      C_item[key] = elnsum(C_item[key], temp_C[key])
  
  return (FW,BW,Gamma, C_item, EC)

################################
def MstepObservation(data_items, ECList, CList, GammaList, Xg, Yg, M_old, O_old):
  """M-step for latent HMM."""
  global FLAG_VB
  global next_states
  global prev_states

  # initialize data structures for observation
  O = {}
  for x in Xg:
    O[x] = (0.0,0.0)
  # initialize data structures for co-occ probs
  M = defaultdict(dict);
  for x in Xg:
    for y in Yg:
      M[x][y] = 0.0;

  # iterate over the data items and estimate aggregated M, C, and O
  for index in range(len(data_items)):
    nouns = data_items[index].nouns
    tags = data_items[index].tags
    Gamma = GammaList[index]
    EC = ECList[index]
    S, SP, EC = IBMProbabilitiesObservation(nouns, tags, M_old, O_old)
    m = len(nouns)
    n = len(tags)
    for i in range(1,m+1):
      for j in range(1,n+1):
        nn = nouns[i-1]
        tt = tags[j-1]
        for key in SP[i][j]:
          prob = SP[i][j][key]/S[i][j]
          noun_set = ListOfNounsFromKey(key)
          for noun in nn:
            if noun in noun_set:
              O[noun] = (O[noun][0], O[noun][1] + prob * Gamma[i][j])
            else:
              O[noun] = (O[noun][0] + prob * Gamma[i][j], O[noun][1])
          for noun in noun_set:
            for tag in tt:
              M[noun][tag] = M[noun][tag] + (Gamma[i][j] * prob * EC[i][j][key][noun][tag])

  # Normalize the O values
  for x in Xg:
    Z = O[x][0] + O[x][1]
    O[x] = (O[x][0]/Z, O[x][1]/Z)
    if O[x][0] < 1e-40:
      O[x]  = (1e-40, O[x][1])
    if O[x][1] < 1e-40:
      O[x]  = (O[x][0], 1e-40)

  # Normalize the M values
  dir_alpha = 1.0
  if FLAG_VB == False:
    dir_alpha = 0.0
  K = float(len(Yg))
  K_alpha = K * dir_alpha
  
  for x in Xg:
    Z = 0.0
    for y in Yg:
      Z = Z + M[x][y];
    
    for y in Yg:
      if FLAG_VB == True:
        if(math.exp(digamma(Z + K_alpha)) < 1e-40):
          M[x][y] = 1e-40
        else:
          M[x][y] = math.exp(digamma( M[x][y]+dir_alpha) )  / math.exp(digamma(Z + K_alpha));
      else:
        M[x][y] = M[x][y]/ Z;
      if M[x][y] < 1e-40:
        M[x][y] = 1e-40
  
  # udpate the C probabilities
  # Normalize the jump probabilities, this should be moved to the M-step part
  Z_c = eln(0.0)
  # initialize to zero
  C = {}
  for i in next_states:
    C[str(i)] = eln(0.0)
  for index in range(len(data_items)):
    C_item = CList[index]
    for i in next_states:
      key = str(i)
      C[key] = elnsum(C[key], C_item[key])
      Z_c = elnsum(Z_c, C_item[key])
  Z_c_inv = Z_c * -1.0
  
  for i in next_states:
    key = str(i)
    C[key] = eexp(elnproduct(C[key], Z_c_inv))
    if (C[key] < 1e-40 ):
      C[key] = 1e-40
  
  print("Current C value:")
  print(C)
  return M, C, O
################################
def ApplyLatentObservationEM(data_items, M, C, O, Xg, Yg, nIterations):
  for iter in range(nIterations):
    ECList = []
    GammaList = []
    CList = []
    print("Iteration: " + str(iter))
    # E-step, by iterating over all the data items
    for item in data_items:
      FW, BW, Gamma, C_item, EC = forward_backward_observation(item.nouns, item.tags, M, C, O)
      GammaList.append(Gamma)
      ECList.append(EC)
      CList.append(C_item)
    # M-step
    M, C, O = MstepObservation(data_items, ECList, CList, GammaList, Xg, Yg, M, O)
  return M, C, O
################################
def JointViterbiLatentObservation(nouns, tags, X, Y, M, C, O):
  """Joint Viterbi Decoding using Latent HMM for single video."""
  m = len(nouns)
  n = len(tags)
  D = [[1e-200 for j in range(n+1)] for i in range(m+1)]
  B= [[(-1,-1) for j in range(n+1)] for i in range(m+1)]
  S, SP, EC = IBMProbabilitiesObservation(nouns, tags, M, O)
  Match = defaultdict(dict)
  
  for x in X:
    for y in Y:
      Match[x][y] = 0.0
  
  # Initialize the viterbi scores
  for i in range(1,m+1):
    D[i][0] = -1e200;
  
  for j in range(1,n+1):
    D[0][j] = -1e200;
  
  
  # Viterbi dynamic programming
  for j in range(1,n+1):
    for i in range(1,m+1):
      nn = nouns[i-1];
      tt = tags[j-1];
      
      BP = GetBackPointers(i,j, m, n)
      max_val = -1e2000
      jump_size = -1000
      for k in range(len(BP)):
        i1 = BP[k][0]
        j1 = BP[k][1]
        
        jump_size = i - i1
        temp_val = D[i1][j1] + math.log(S[i][j]+1e-200) + math.log(C[str(jump_size)]+1e-200)
        if  temp_val >= max_val:
          max_val = temp_val
          B[i][j] = (i1,j1)
          jump_size = i-i1
      D[i][j]  = max_val
  return D,B, Match

############################
def ViterbiAlignmentAccuraciesLatentHMM(data_items, M, C, O):
  """Viterbi decoding for all the data items using Latent HMM."""
  accuracies = []
  avg_acc = 0.0
  nItems = len(data_items)
  for item in data_items:
    D, B, Match = JointViterbiLatentObservation(item.nouns, item.tags, item.X, item.Y, M, C, O)
    bestAlignment = GetBestAlignment(item.nouns, item.tags, B, D)
    acc = StepwiseAlignmentAccuracy(bestAlignment, item.step_index_avl, item.step_index_nouns)
    print("Accuracy: %f" % acc)
    accuracies.append(acc)
    avg_acc = avg_acc + acc/float(nItems)
  print "Avg accuracy HMM: %f" % avg_acc
  return accuracies, avg_acc