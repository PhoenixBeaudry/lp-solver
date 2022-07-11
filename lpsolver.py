import numpy
import sys
debug = False


# Read in LP from stdin
lp = []
for line in sys.stdin:
    coeff_line = []
    for coeff in line.split():
        coeff_line.append(float(coeff))
    lp.append(coeff_line)

if (debug):
    print(lp)


