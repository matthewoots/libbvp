#!/usr/bin/env python3

import math
import sys 

# Delta wing example for ZoHD dart 250g 
# M = 0.2565, T = 5 in 4s = 0.8s, b = 0.325, L = 0.35
# I = 0.012311

# https://physics.stackexchange.com/questions/147452/calculating-the-moment-of-inertia-in-bifilar-pendulums

# M = mass
# T = rotational time period for one rotating oscillation 
# b = length between where string attaches to rod and centre of gravity
# L = length of string that suspends rod

if len(sys.argv) < 4:
    # print('error %s is %d years old' % (name, age))
    print('error sys.argv < 4')
    sys.exit()

M = float(sys.argv[1])
T = float(sys.argv[2])
b = float(sys.argv[3])
L = float(sys.argv[4])
g = 9.81

I = (M * g * math.pow(T,2) * math.pow(b,2)) / (4 * math.pow(3.1415,2) * L)

print('Moments of inertia for axis = %lf ' % (I))