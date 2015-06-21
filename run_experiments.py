#!/usr/bin/python

"""
  To run all the experiments using this script, use the command:
    > python run_experiments.py <path_to_output_directory>
"""

import sys
import math
import operator
import random
import os


def GetCommandLine(file, hmm_flag, anvil_flag):
  """Construct the command line argument to run the alignment code."""
  output_file = "./%s/Results_%s_%s_%s.txt" % (sys.argv[1],file, hmm_flag, anvil_flag)
  exec_line = "python %s.py %s %s > %s" % (file, hmm_flag, anvil_flag, output_file)
  return exec_line

############################
def main():
  code_list = ["hmm_multivideos" ,"latent_hmm_multivideos"]
  hmm_flag_list = ["True", "False"]
  anvil_flag_list = ["True", "False"]

  for code in code_list:
    for hmm_flag in hmm_flag_list:
      for anvil_flag in anvil_flag_list:
        exec_line = GetCommandLine(code, hmm_flag, anvil_flag)
        print exec_line
        os.system(exec_line)

if __name__ == "__main__":
  main()
