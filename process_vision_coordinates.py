#!/usr/bin/python
#
# Author: Iftekhar Naim, CS Dept., University of Rochester.
# Email: inaim@cs.rochester.edu
#

"""
  This script converts the raw vision tracks to Anvil format.
  It is not directly related to the main HMM alignment code.
  If you have just started looking at the code, you can skip it.
"""
import sys
import math
import operator
import random

from os import listdir
from os.path import isfile, join
import matplotlib
import pylab


# load custom class
from timing_info import TimeInterval
from input_file_processor import *

special_files = {}

global vision_dir_path;

# get the additional special data
def GetSpecialCoordinates(video_index):
  in_hand_list = ()
  unique_objects = ()
  if (video_index == "cell06") | (video_index == "cell08") |  (video_index == "cell21"):
    (in_hand_list, unique_objects) = GetSpecialInHandLists(vision_dir_path + "post_processed_files/obj_hands2_kalman.txt")
    
    (in_hand_list2, unique_objects2) = GetSpecialInHandLists(vision_dir_path + "post_processed_files/obj_yellow2_kalman.txt")
    
    for ii in range(len(in_hand_list2)):
      in_hand_list.append(in_hand_list2[ii])
      unique_objects.append(unique_objects2[ii])
  
  
  elif (video_index == "llgm10") | (video_index == "llgm14"):
    (in_hand_list, unique_objects) = GetSpecialInHandLists(vision_dir_path + "post_processed_files/obj_hands2_kalman.txt")
  
  #        (in_hand_list2, unique_objects2) = GetSpecialInHandLists(vision_dir_path + "post_processed_files/obj_yellow2_kalman.txt")
  #
  #        for ii in range(len(in_hand_list2)):
  #            in_hand_list.append(in_hand_list2[ii])
  #            unique_objects.append(unique_objects2[ii])
  
  
  elif (video_index == "ypad11") | (video_index == "ypad15"):
    (in_hand_list, unique_objects) = GetSpecialInHandLists(vision_dir_path + "post_processed_files/obj_hands2_kalman.txt")
  
  
  print(len(in_hand_list))
  print(unique_objects)
  
  
  return(in_hand_list, unique_objects)

def GetSpecialInHandLists(filepath):
  unique_objects = []
  f = open(filepath, 'r')
  lines = f.readlines()
  f.close()
  nFrames = len(lines)
  
  print(filepath + " : " + str(nFrames) )
  
  for i in range(len(lines)):
    lines[i] = lines[i].lstrip().rstrip()
    parts = lines[i].split(" ")
    if parts[-1] == "UNKNOWN":
      continue
    if parts[-1] not in unique_objects:
      unique_objects.append(parts[-1])
  in_hand_list = []
  
  print(unique_objects)
  raw_input("")
  
  for obj in unique_objects:
    in_hand = [0] * nFrames
    for i in range(len(lines)):
      parts = lines[i].split(" ")
      if parts[-1] == obj:
        in_hand[i] = 1
      else:
        in_hand[i] = 0
    
    in_hand_list.append(in_hand)
  
  return (in_hand_list, unique_objects)

############################
def GetCoordinateList(vision_dir_path):
  coordinate_list = []
  
  tempfiles = [ f for f in listdir(vision_dir_path) if isfile(join(vision_dir_path,f)) ]
  
  files = []
  for f in tempfiles:
    # ignore any file that starts with a "."
    if f[0] == '.':
      continue
    files.append(f)
  
  
  print(files)
  
  
  for filename in files:
    f = open(vision_dir_path + filename, 'r')
    lines = f.readlines()
    f.close()
    
    coordinates = Get3DCoordinates(lines)
    print(filename)
    print(len(coordinates))
    
    coordinate_list.append(coordinates)
    
    raw_input()
  
  nFrames = len(coordinate_list[0])
  nObjects = len(coordinate_list)
  
  print("nFrames = " + str(nFrames))
  print("nObjects = " + str(nObjects))
  return (coordinate_list, nFrames, nObjects, files)

############################
def Get3DCoordinates(lines):
  coordinates = []
  for i in range(len(lines)):
    parts = lines[i].split(" ")
    coordinates.append( (float(parts[0]),float(parts[1]),float(parts[2])) )
  return coordinates

def L2Distance(a,b):
  if not (len(a) == len(b)):
    return -1.0
  
  dist = 0.0
  for i in range(len(a)):
    dist = dist + (a[i] - b[i])**2
  
  return dist

############################
def GetFrameTimings(video_index, nFrames):
  timings = []
  
  path = "./frames/frames_" + video_index + ".txt";
  
  f = open(path, "r")
  lines = f.readlines()
  
  count = -1
  
  start_time = 0.0
  for line in lines:
    count = count + 1
    #ni3-000001-138186986811-rgb.png
    parts = line.lstrip().rstrip().split("-")
    msec = int(parts[2])
    if count == 0:
      start_time = msec
    
    
    if (video_index == "cell21"):
      if not((count%5)==0):
        continue
      divisor = 1000.0
    else:
      divisor = 100.0
    
    if  (video_index == "cell22") or (video_index == "cell23") or (video_index == "cell24") or (video_index == "cell24") or (video_index == "cell30"):
      divisor = 1000.0

    timeval = float(msec - start_time) / divisor
    
    print(msec, timeval, count, path)
    
    timings.append( timeval )
    
    if count == nFrames-1:
      break
  
  return timings

def GetTimeIntervalWithStepGroundTruth(video_index):
  true_anvil_path = './anvil_tags/' + video_index.upper() + '.anvil.txt';
  
  
  f = open(true_anvil_path, "r");
  
  lines = f.readlines()
  f.close()
  
  start_time = 0.0
  
  prev_step = 1
  
  ground_truth_timing = []
  
  for i in range(1, len(lines)):
    parts = lines[i].split(" ")
    
    if len(parts) < 4:
      continue
    last_part = parts[3]
    last_parts = last_part.split(",")
    if len(last_parts) < 2:
      continue
    
    
    new_step = int(last_parts[-1])
    if new_step == prev_step:
      continue
    else:
      end_time = float(parts[1])
      ground_truth_timing.append((start_time, end_time,prev_step))
      start_time = end_time
      prev_step = new_step
  
  ground_truth_timing.append((start_time, 1000000, new_step))
  
  
  return ground_truth_timing


def GetStep(start_time, end_time, ground_truth_timing):
  for i in range(len(ground_truth_timing)):
    if (start_time >= ground_truth_timing[i][0]) & (start_time <= ground_truth_timing[i][1]):
      return ground_truth_timing[i][2]
############################
def main():
  global vision_dir_path;
  vision_dir_path = sys.argv[1]
  # if path doesn't end with a "/", then add
  if not(vision_dir_path[-1] == '/'):
    vision_dir_path = vision_dir_path + "/";
  
  path_parts = vision_dir_path.split("/")
  video_index = path_parts[-2]
  
  print(video_index)

  # Threshold is 0.02 for older data, and 0.01 for the new data
  th = 0.004
  print(vision_dir_path)
  
  (coordinate_list, nFrames, nObjects, files) = GetCoordinateList(vision_dir_path)
  

  timings = GetFrameTimings(video_index, nFrames)
  print(timings)

  ground_truth_timing = GetTimeIntervalWithStepGroundTruth(video_index)
  print(ground_truth_timing)
  
  in_hand_list = []
  
  for i in range(2, nObjects):
    
    distance_left = [0.0] * nFrames
    distance_right = [0.0] * nFrames
    time = [0.0] * nFrames
    
    in_hand = [0] * nFrames
    
    for j in range(0,nFrames):
      distance_left[j] = L2Distance(coordinate_list[0][j], coordinate_list[i][j])
      
      distance_right[j] = L2Distance(coordinate_list[0][j], coordinate_list[i][j])
      time[j] = float(j) / 6.0
      
      if (distance_left[j] < th) | (distance_right[j] < th):
        in_hand[j] = 1
    
    
    
    fig = matplotlib.pyplot.figure()
    matplotlib.pyplot.scatter(time,distance_left)
    matplotlib.pyplot.ylim(0, 0.4)
    fig.suptitle('Left hand ' + files[i] , fontsize=20)
    matplotlib.pyplot.show()
    
    
    fig = matplotlib.pyplot.figure()
    matplotlib.pyplot.scatter(time,distance_right)
    matplotlib.pyplot.ylim(0, 0.4)
    fig.suptitle('Right hand ' + files[i] , fontsize=20)
    matplotlib.pyplot.show()
    
    
    #print(in_hand)
    
    in_hand_list.append(in_hand)
  
  print("Initial number of items:" + str(len(in_hand_list)))
  # get few more additional in_hand lists
  (in_hand_list_2, files2) = GetSpecialCoordinates(video_index)
  
  for i in range(len(files2)):
    files.append(files2[i])
    in_hand_list.append(in_hand_list_2[i])
  
  fout = open(vision_dir_path + "../" + video_index + ".vision.anvil.txt", "w")
  fout.write("ID start end token\n");
  # iterate th  in_hand_list
  count = 0
  
  prev_objects = []
  start_time = -1
  end_time = 0.0
  for t in range(nFrames):
    
    objects = []
    for j in range(len(in_hand_list)):
      if in_hand_list[j][t] == 1:
        objects.append(files[j+2])
    
    if not(prev_objects == objects) | (t == nFrames - 1):
      #            print("----")
      #            print(prev_objects)
      #            print(objects)
      #            raw_input(" ")
      
      # for the first time found any object touching hands
      if start_time == -1:
        start_time = timings[t]
        prev_objects = objects
        continue
      
      # A hard threshold to ignore spurious events having less than 0.2 sec duration
      if (timings[t] - start_time) < 0.2:
        continue
      # write down the events
      if(len(prev_objects) > 0):
        fout.write("f0t0e0 " + str(start_time) + " " + str(timings[t]) + " ");
        for ll in range(len(prev_objects)):
          fout.write(prev_objects[ll]+",")
        
        step = GetStep(start_time, timings[t], ground_truth_timing)
        fout.write(str(step) + "\n")
      
      start_time = timings[t]
      prev_objects = objects
  
  fout.close()
  print(vision_dir_path + "../" + video_index + "_anvil.txt")
#        raw_input("")
if __name__ == "__main__":
  main()
