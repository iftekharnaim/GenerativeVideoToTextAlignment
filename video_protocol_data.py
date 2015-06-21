#!/usr/bin/python
#
# Author: Iftekhar Naim, CS Dept., University of Rochester.
# Email: inaim@cs.rochester.edu
#
"""A class for storing Video and Protocol Data."""

class VideoProtocolData:
  """The class containing data for a (video, protocol) pair."""
  
  def __init__(self, nouns, tags, verbs, step_index_nouns, step_index_avl, lines_text, lines_avl, X, Y, Z, timeX, timeY, timing_info_anvil):
    self.nouns = nouns
    self.tags = tags
    self.verbs = verbs

    self.m = len(nouns)
    self.n = len(tags)
    
    self.step_index_nouns = step_index_nouns
    self.step_index_avl = step_index_avl
    self.lines_text = lines_text
    self.lines_avl = lines_avl

    self.X = X
    self.Y = Y
    self.Z = Z
    
    self.timeX = timeX
    self.timeY = timeY
    self.timing_info_anvil = timing_info_anvil
    self.video_timing_gaps = self.ReadVideoChunkTimeGaps(timing_info_anvil)

  ############################
  def ReadVideoChunkTimeGaps(self, timing_info_anvil):
    nChunks = len(timing_info_anvil)
    time_gaps = []
    time_gaps.append(0.0)
    for i in range(1, nChunks):
      prev_tinfo = timing_info_anvil[i-1]
      tinfo = timing_info_anvil[i]
      interval = tinfo.start - prev_tinfo.end
      time_gaps.append(interval)
    assert len(time_gaps) == nChunks
    return time_gaps
