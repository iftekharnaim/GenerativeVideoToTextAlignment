#!/usr/bin/python
#
# Author: Iftekhar Naim, CS Dept., University of Rochester.
# Email: inaim@cs.rochester.edu
#

"""A file containing all the global variables."""

# Whether HMM2 (more flexible transitions). Default to monotonic.
FLAG_HMM2 = False
# Whether ANVIL data (true) or vision data (false)
FLAG_ANVIL = True
# Whether to perform variational bayes. Usually set to False.
FLAG_VB = False
# next and prev states
next_states = range(0,2,1)
prev_states = range(-1,1,1)
# Total number of EM iterations
nIterationsEM = 50

# datasets
# define a list of protocols and corresponding videos
protocols_anvil = ["CELL", "CELL", "LLGM", "LLGM", "YPAD", "YPAD"]
anvil_tags = ["CELL06.anvil.txt", "CELL21.anvil.txt",
              "LLGM10.anvil.txt", "LLGM14.anvil.txt",
              "YPAD11.anvil.txt", "YPAD15.anvil.txt"]

protocols_vision = ["CELL", "CELL", "LLGM", "LLGM", "YPAD", "YPAD"]
vision_tags = [
   "cell06.vision.anvil.txt", "cell21.vision.anvil.txt",
   "llgm10.vision.anvil.txt", "llgm14.vision.anvil.txt",
   "ypad11.vision.anvil.txt", "ypad15.vision.anvil.txt"]