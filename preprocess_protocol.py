"""A simple preprocessing code for protocols."""
# In order to parse, the BLLIP parser code requires
# the data in a specific format.
# E.g., <s> at the beginning and </s> at the end of each line.
# This simple script does this preprocessing.

import sys
# get the input file name
filename = sys.argv[1];
# open the file
f = open(filename)
lines = f.readlines()
f.close()

# open a new preprocessed file
outputfile = "preprocessed_file" + filename;
fw = open(outputfile, 'w');

# print the files
for line in lines:
  line = line.lstrip()
  if (len(line)== 0):
    continue
  print("<s> " + line.rstrip() + " </s>")