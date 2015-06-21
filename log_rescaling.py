#!/usr/bin/python
#
# Author: Iftekhar Naim, CS Dept., University of Rochester.
# Email: inaim@cs.rochester.edu
#

import sys
from collections import defaultdict
import math
import operator

# Check if two floating point numbers are too close
def isnear(x,y):
    if math.fabs(x-y) < 1e-100:
        return True

# The extended exponential function
# Handles the NaN input
def eexp(x):
    if math.isnan(x):
        return 0.0;
    return math.exp(x);


# The extended logarithm function
# Handles the NaN input
def eln(x):
    if isnear(x, 0.0):
        return float('nan');
    
    if x < 0.0:
        print('Error. Negative value inside logarithm:' + str(x))

    return math.log(x);


# Compute sums in log domain.
# Inputs: x = log(a), and y = log(b)
# Output: log(a+b)
def elnsum(x,y):
    if (math.isnan(x)) | (math.isnan(y)):
        if math.isnan(x):
            return y
        elif math.isnan(y):
            return x
    else:
        if (x > y):
            return (x + eln(1.0 + math.exp(y-x)))
        else:
            return (y + eln(1.0 + math.exp(x-y)))
    
# Compute the product in log domain
def elnproduct(x,y):
    if (math.isnan(x)) | (math.isnan(y)):
        return float('nan')
    else:
        return (x+y)

############################################
# For unit testing
def main():
    
    print(eln(math.exp(10)));
    print(eln(0.0));

    a = eln(0.0)
    b = eexp(a)
    print(b)
    print(eexp(eln(10)))
                    
    sum = elnsum(eln(50) , eln(50))
    print(sum)
    print(eln(100.0))
                    
    prod = elnproduct(eln(50) , eln(50))
    print(prod)
    print(eln(2500))

if __name__ == "__main__":
    main()




