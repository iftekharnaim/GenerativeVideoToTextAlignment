class TimeInterval:
    """The class containing timing info for anvil tags and nouns"""

 
    def __init__(self, start = -1.0, end = 1.0, mean = -1.0):
        self.start = start # The starting time of the interval
        self.end = end # The ending of the time interval
        self.mean = mean # The center of the interval
        
    def PrintTimeInterval(self):
        print('start: ' + str(self.start) + "  end:" + str(self.end) + "  mean:" + str(self.mean));
